import re
from Bio import AlignIO, SeqIO
from collections import defaultdict
from opsintools.scripts import hmmer

class Hmmer:

    @staticmethod
    def parse_hmm_name_len(line):
        """Parse the profile length line."""
        hmm_name, hmm_len = re.findall(r'([^ ]+) +\[M=(\d+)\]', line)[0]
        return hmm_name, int(hmm_len)

    @staticmethod
    def parse_complete_data(fh):
        """Parse the full scores block.

        The block ends with an empty line.
        """
        inclusion = True
        values = {}
        while (line := next(fh).strip()) != "":
            if 'inclusion threshold' in line:
                inclusion = False
            if not line.startswith('--') and not line.startswith('E-value'):
                full_E_value, full_score, full_bias, best_E_value, best_score, best_bias, exp, N, Sequence, *rest = line.split()
                values[Sequence] = {
                    "full": {
                        "evalue": float(full_E_value),
                        "score": float(full_score),
                        "bias": float(full_bias)
                    },
                    "best": {
                        "evalue": float(best_E_value),
                        "score": float(best_score),
                        "bias": float(best_bias)
                    },
                    "exp": float(exp),
                    "num_domains": int(N),
                    "included": inclusion
                }
        return values

    @staticmethod
    def parse_seq_name_desc(line):
        """Parse sequence name and description."""
        seq_name, seq_desc = re.findall(r'>>\s*(\S+)\s*(.*)', line)[0]
        return seq_name, seq_desc.strip()

    @staticmethod
    def completeness(compl):
        return [ compl[0] == '[', compl[1] == ']' ]

    @staticmethod
    def parse_domain_data(fh):
        """Parse a block for one input sequence.

        The block starts with ">>" and
        ends with last domain's last alignment block.
        """
        domains = []
        no_domains = False

        # The scores section starts immediately and ends with an empty line
        while (line := next(fh).strip()) != '' and not no_domains:
            no_domains = line.startswith("[No individual domains")
            if not line.startswith('#') and not line.startswith('--') and not no_domains:
                num, passed, score, bias, c_evalue, i_evalue, hmm_from, hmm_to, hmm_compl, ali_from, ali_to, ali_compl, env_from, env_to, env_compl, acc = line.split()
                vals = {
                    "num": int(num),
                    "passed": passed == "!",
                    "score": float(score),
                    "bias": float(bias),
                    "c_evalue": float(c_evalue),
                    "i_evalue": float(i_evalue),
                    "hmm": {
                        "from": int(hmm_from),
                        "to":   int(hmm_to),
                        "complete": Hmmer.completeness(hmm_compl),
                        "seq": ""
                    },
                    "ali": {
                        "from": int(ali_from),
                        "to":   int(ali_to),
                        "complete": Hmmer.completeness(ali_compl),
                        "seq": ""
                    },
                    "env": {
                        "from": int(env_from),
                        "to":   int(env_to),
                        "complete": Hmmer.completeness(env_compl)
                    },
                    "acc": float(acc),
                    "RF": "",
                    "PP": "",
                    "matches": ""
                }
                domains.append(vals)

        alignments_for_each_domain = next(fh) # skip this line

        i = -1
        while i < len(domains) - 1:
            line = next(fh).strip()
            if line.startswith('=='):
                i += 1
                hmm_right = -1
                while hmm_right < domains[i]['hmm']['to']:
                    block = Hmmer.parse_aln_block(fh)
                    domains[i]["RF"] += block['RF']
                    domains[i]["PP"] += block['PP']
                    domains[i]["hmm"]["seq"] += block['hmm_seq']
                    domains[i]["ali"]["seq"] += block['ali_seq']
                    domains[i]["matches"] += block['matches']
                    if block['hmm_right'] != "-":
                        hmm_right = int(block['hmm_right'])
        return domains

    @staticmethod
    def parse_aln_block(fh):
        """Parse profile-query alignment (sub)block.

        The block ends with an empty line.
        """
        vals = {}

        line = next(fh)
        fields = re.split(' +', line.strip(), 1)
        if fields[1] == 'RF':
            vals['RF'] = fields[0]
            line = next(fh)
        else:
            vals['RF'] = ""

        fields = re.split(" +", line.strip())
        vals['hmm_left']  = fields[-3]
        vals['hmm_seq']   = fields[-2]
        vals['hmm_right'] = fields[-1]

        line = next(fh)
        vals['matches'] = line.rstrip('\n\r')[-len(vals['hmm_seq']):]

        line = next(fh)
        fields = re.split(" +", line.strip())
        vals['ali_left']  = fields[-3]
        vals['ali_seq']   = fields[-2]
        vals['ali_right'] = fields[-1]

        line = next(fh)
        fields = re.split(' +', line.strip(), 2)
        vals['PP'] = fields[0]

        line = next(fh).strip()
        assert line == '', f"An empty line expected, got {line}"

        return vals

    def __init__(self, profile_file, profile_cons, search = None, align = None):
        self.matches = []
        profile_file = str(profile_file)
        if search:
            self.parse_hmmsearch(search, profile_cons, profile_file)
        if align:
            self.parse_hmmalign(align, profile_cons, profile_file)

    def parse_hmmsearch(self, hmmsearch_file, profile_cons, profile_file):
        hmm_name = hmm_len = None
        complete_data = {}
        if hmmsearch_file:
            with open(hmmsearch_file) as file:
                for line in file:
                    if line.startswith('Query:'):
                        hmm_name, hmm_len = Hmmer.parse_hmm_name_len(line)
                    elif line.startswith('Scores'):
                        complete_data = Hmmer.parse_complete_data(file)
                    elif line.startswith('>>'):
                        seq_name, seq_desc = Hmmer.parse_seq_name_desc(line)
                        if seq_name not in complete_data:
                            raise ValueError(f"Sequence does not appear in the complete sequences table: {seq_name}")
                        match = complete_data[seq_name]
                        match["description"] = seq_desc
                        match["domains"] = Hmmer.parse_domain_data(file)
                        match["seq_name"] = seq_name
                        match["hmm_name"] = hmm_name
                        match["hmm_len"] = hmm_len
                        match["program"] = "hmmsearch"
                        match["file_name"] = str(hmmsearch_file)
                        match["profile_file"] = profile_file
                        match["profile_cons"] = profile_cons
                        self.matches.append(match)

    @staticmethod
    def parse_sto_record(record):
        ali_aln = []
        pp_aln = []
        ali_pos = 0
        hmm_pos = 0
        ali_from = ali_to = None
        hmm_from = hmm_to = None
        for aa, pp in zip(record.seq, record.letter_annotations['posterior_probability']):
            if aa != '.' and pp != '.': # not a space-filler
                is_aligned_residue = aa.isupper()
                is_aligned = is_aligned_residue or aa == '-'
                is_residue = aa != '.'
                if is_aligned:
                    hmm_pos += 1
                    if hmm_from is None:
                        hmm_from = hmm_pos
                if is_residue:
                    ali_pos += 1
                    if is_aligned_residue:
                        if ali_from is None:
                            ali_from = ali_pos
                        if hmm_from is None:
                            hmm_from = hmm_pos
                        ali_to = ali_pos
                        hmm_to = hmm_pos
                ali_aln.append(aa)
                pp_aln.append(pp)
        hmm_len = hmm_pos
        return ''.join(ali_aln), ''.join(pp_aln), ali_from, ali_to, hmm_from, hmm_to, hmm_len

    def parse_hmmalign(self, sto_file, profile_cons, profile_file):
        align = AlignIO.read(sto_file, "stockholm")
        for record in align:
            ali_seq, pp_seq, ali_from, ali_to, hmm_from, hmm_to, hmm_len = Hmmer.parse_sto_record(record)
            match = {
                "seq_name": record.id,
                "hmm_name": "",
                "hmm_len": hmm_len,
                "description": record.description,
                "full": { "evalue": -1, "score": -1, "bias": -1 },
                "best": { "evalue": -1, "score": -1, "bias": -1 },
                "exp": -1,
                "num_domains": 1,
                "included": True
            }
            match["domains"] = [ {
                "num": 1,
                "passed": None, "score": None, "bias": None, "c_evalue": None, "i_evalue": None,
                "hmm": {
                    "seq": "",
                    "from": hmm_from,
                    "to": hmm_to
                },
                "ali": {
                    "seq": ali_seq,
                    "from": ali_from,
                    "to": ali_to
                },
                "PP": pp_seq
            } ]
            match["program"] = "hmmalign"
            match["file_name"] = str(sto_file)
            match["profile_file"] = profile_file
            match["profile_cons"] = profile_cons
            self.matches.append(match)

    @staticmethod
    def decode_prob(pp_res):
        if pp_res == "*":
            return 10
        if pp_res.isdigit():
            return int(pp_res)
        if pp_res == ".":
            return -1
        raise ValueError(f"Unrecognized posterior probability {pp_res}")

    @staticmethod
    def encode_prob(pp_val):
        if pp_val == 10:
            return "*"
        if -1 < pp_val < 10:
            return str(pp_val)
        if pp_val == -1:
            return "."
        raise ValueError(f"Unrecognized posterior probability value {pp_val}")

    def chain_domains(self, records, max_gap = 100):
        for match in self.matches:
            seq_name = match["seq_name"]
            if seq_name not in records:
                raise ValueError(f"{seq_name} not in fasta")
            record = records[seq_name]
            profile_cons = match["profile_cons"]
            domains = match["domains"]
            merged = hmmer.chain_local_alignments(str(record.seq), profile_cons, domains, max_gap = max_gap)
            match['domains'] = merged

class HmmerContainer(Hmmer):
    def __init__(self):
        self.matches = []

    def resolve(self):
        best_matches = {}
        for match in self.matches:
            seq_name = match["seq_name"]
            if seq_name not in best_matches or best_matches[seq_name]["full"]["score"] < match["full"]["score"]:
                best_matches[seq_name] = match
        self.matches = list(best_matches.values())
        return self

    def __add__(self, hmmer2: Hmmer):
        self.matches += hmmer2.matches
        return self
