#@title Output

import json
from IPython.display import display, HTML

# --- Visualization ---
MAX_ALNS = 10
CHUNK_SIZE = 100 # Alignment block width

zappo = {
    'I': '#FFAFAF', 'L': '#FFAFAF', 'V': '#FFAFAF', 'A': '#FFAFAF', 'M': '#FFAFAF',
    'F': '#FFC800', 'W': '#FFC800', 'Y': '#FFC800',
    'K': '#6464FF', 'R': '#6464FF', 'H': '#6464FF',
    'D': '#FF0000', 'E': '#FF0000',
    'S': '#00FF00', 'T': '#00FF00', 'N': '#00FF00', 'Q': '#00FF00',
    'P': '#FF00FF', 'G': '#FF00FF',
    'C': '#FFFF00'
 }

def get_color_rgba(aa, score):
    hex_color = zappo.get(aa.upper(), '#FFFFFF')
    alpha = score / 10 if score is not None else 0.1
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)
    return f"rgba({r}, {g}, {b}, {alpha})"

def build_number_track(sequence, start_idx):
    track = [' '] * len(sequence)
    current_idx = start_idx
    for j, char in enumerate(sequence):
        if char.isalpha():
            current_idx += 1
            if current_idx % 10 == 0:
                num_str = str(current_idx)
                for k, digit in enumerate(num_str):
                    if j + k < len(track):
                        track[j + k] = digit
                    else:
                        track.append(digit)

    return "".join(track), current_idx

def build_tm_track(sequence, start_idx, tms):
    tm_at_pos = [None] * len(sequence)
    current_idx = start_idx

    for j, char in enumerate(sequence):
        is_alpha = char.isalpha()
        if is_alpha:
            current_idx += 1

        for tm_name, coords in tms.items():
            if is_alpha:
                if coords[0] <= current_idx <= coords[1]:
                    tm_at_pos[j] = tm_name
                    break
            else:
                if coords[0] <= current_idx < coords[1]:
                    tm_at_pos[j] = tm_name
                    break

    tm_html = []
    j = 0
    while j < len(sequence):
        tm = tm_at_pos[j]
        if tm is None:
            tm_html.append(" ")
            j += 1
        else:
            k = j
            while k < len(sequence) and tm_at_pos[k] == tm:
                k += 1
            block_len = k - j

            label = tm if block_len >= len(tm) else tm[:block_len]
            text_content = label.ljust(block_len, " ")

            tm_html.append(f'<span style="background-color: steelblue; color: white; font-weight: bold;">{text_content}</span>')
            j = k

    return "".join(tm_html), current_idx

def parse_score(score_str):
    if score_str == "*":
        return 10
    if score_str.isdigit():
        return int(score_str)
    return None

def generate_html_report(output_path, data_path, def_profile = None, num_alns = MAX_ALNS):
    json_paths = list(output_path.glob("*/opsinmap.json"))[:num_alns]
    output_html_path = output_path / "alignment.html"

    html_out = [
        "<html><head><style>",
        "body { font-family: monospace; font-size: 14px; line-height: 0.85; }",
        ".chunk { margin-bottom: 24px; white-space: pre; }",
        ".header { font-weight: bold; font-size: 16px; margin-top: 15px; margin-bottom: 10px; border-bottom: 1px solid #ccc; padding-bottom: 5px; }",
        ".row-label { display: inline-block; width: 60px; font-weight: bold; color: #555; }",
        ".match-char { color: #666; font-weight: bold; }",
        ".scrollable-wrapper { max-height: 500px; overflow-y: auto; padding-right: 10px; }",
        "</style></head><body>",
        '<div class="scrollable-wrapper">'
    ]

    ref_data = {}
    for ref_path in data_path.glob('*/ref.json'):
        profile = ref_path.parent.name
        ref_data[profile] = json.loads(ref_path.read_text())

    for json_path in json_paths:
        query_data = json.loads(json_path.read_text())

        for domain in query_data:
            alignment = domain["alignment"]
            profile = domain.get("profile", def_profile)
            assert profile in ref_data, f"{profile} not found in {data_path}"
            tms = ref_data[profile]["tms"]
            q_seq_raw = alignment["query"]
            r_seq_raw = alignment["ref"]
            q_score_raw = alignment["query_score"]
            r_score_raw = alignment["ref_score"]

            q_seq = []
            r_seq = []
            q_score = []
            r_score = []

            for qa, ra, qs, rs in zip(q_seq_raw, r_seq_raw, q_score_raw, r_score_raw):
                if qa != '-' or ra != '-':
                    q_seq.append(qa)
                    r_seq.append(ra)
                    q_score.append(parse_score(qs))
                    r_score.append(parse_score(rs))

            html_out.append(f'<div class="header">Query: {domain["query"]}, reference: {domain["ref"]}, dataset: {profile}</div>')

            q_idx = 0
            r_idx = 0

            for i in range(0, len(q_seq), CHUNK_SIZE):
                q_chunk = q_seq[i : i + CHUNK_SIZE]
                r_chunk = r_seq[i : i + CHUNK_SIZE]
                qs_chunk = q_score[i : i + CHUNK_SIZE]
                rs_chunk = r_score[i : i + CHUNK_SIZE]

                # 1. Query Numbers
                q_num_track, next_q_idx = build_number_track(sequence = q_chunk, start_idx = q_idx)
                html_out.append(f'<div class="chunk"><span class="row-label"></span>{q_num_track}<br>')

                # 2. Query Sequence
                q_html = ['<span class="row-label">Query</span>']
                for j, qa in enumerate(q_chunk):
                    rgba = get_color_rgba(aa = qa, score = qs_chunk[j])
                    q_html.append(f'<span style="background-color: {rgba};">{qa}</span>')
                html_out.append("".join(q_html) + "<br>")

                # 3. Match row
                m_html = ['<span class="row-label"></span>']
                for j in range(len(q_chunk)):
                    qa = q_chunk[j]
                    ra = r_chunk[j]
                    qs = qs_chunk[j]
                    rs = rs_chunk[j]

                    # Check that both scores exist and both chars are alpha
                    if qa.isalpha() and ra.isalpha() and qs is not None and rs is not None:
                        if qa.upper() == ra.upper():
                            m_html.append(f'<span class="match-char">{qa.upper()}</span>')
                        elif zappo.get(qa.upper()) == zappo.get(ra.upper()) and zappo.get(qa.upper()) is not None:
                            m_html.append('<span class="match-char">+</span>')
                        else:
                            m_html.append(" ")
                    else:
                        m_html.append(" ")
                html_out.append("".join(m_html) + "<br>")

                # 4. Reference Sequence
                r_html = ['<span class="row-label">Ref</span>']
                for j, ra in enumerate(r_chunk):
                    rgba = get_color_rgba(aa = ra, score = rs_chunk[j])
                    r_html.append(f'<span style="background-color: {rgba};">{ra}</span>')
                html_out.append("".join(r_html) + "<br>")

                # 5. Reference Numbers
                r_num_track, next_r_idx = build_number_track(sequence = r_chunk, start_idx = r_idx)
                html_out.append(f'<span class="row-label"></span>{r_num_track}<br>')

                # 6. TM Track
                tm_html, _ = build_tm_track(sequence = r_chunk, start_idx = r_idx, tms = tms)
                html_out.append(f'<span class="row-label"></span>{tm_html}</div>')

                q_idx = next_q_idx
                r_idx = next_r_idx

    html_out.append("</div></body></html>")

    html_string = "\n".join(html_out)

    with open(file = output_html_path, mode = "w") as f:
        f.write(html_string)

    display(HTML(data = html_string))
