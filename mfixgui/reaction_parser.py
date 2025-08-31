#!/usr/bin/env python

#from pprint import pprint

class ParseError(Exception):
    pass

# Turn values into numbers, where possible
def to_float(val):
    vtmp = val.lower().replace('d', 'e')
    try:
        return float(vtmp)
    except ValueError:
        return val

class ReactionParser:
    def __init__(self):
        self.state = 0
        self.reactions = {}


    def parse(self, data):
        for line in data.splitlines():
            self.parse_line(line)

    def parse_line(self, line):
        line = line.strip()
        if '!' in line: # comment
            line = line.split('!')[0].strip()
        if not line: #blank
            return
        #print("PARSE STATE=", self.state, line)

        if self.state == 0: # Beginning of reaction
            if '{' not in line:
                raise ParseError('Expected "{" character: %s' % line)
            self.reaction_id = line.split('{',1)[0].strip()
            if not self.valid_reaction_id(self.reaction_id):
                raise ValueError("Invalid reaction name %s" % self.reaction_id)
            if self.reaction_id in self.reactions:
                raise ValueError("Duplicate reaction %s" % self.reaction_id)
            self.reaction = self.reactions[self.reaction_id] = {}
            if line.endswith('}'): # Single-line format
                #Reaction { key = val }
                line = line.split('{', 1)[1].split('}',1)[0]
                key, val = line.split('=',1)
                key = key.strip().lower()
                val = val.strip()
                if key.lower() == 'chem_eq':
                    self.reaction['chem_eq'] = val
                    self.finish_reaction()
                else:
                    raise ParseError("Expected chem_eq: %s" % line)
                self.state = 0
                return
            if not line.endswith('{'):
                raise ParseError('Unexpected text after "{" character: %s'  %line)
            self.state = 1
        elif self.state == 1: # Start of line in reaction body
            if line == '}':
                self.finish_reaction()
                self.state = 0
                return
            if not '=' in line:
                raise ParseError('Expected "=" character: %s' % line)
            self.key = line.split('=',1)[0].strip().lower()
            self.val = line.split('=',1)[1].strip()
            if self.val.endswith(('&','}')):
                self.val = self.val[:-1].strip()
            try:
                self.reaction[self.key] = to_float(self.val)
            except ValueError:
                self.reaction[self.key] = self.val
            if line.endswith('&'):
                self.state = 2
            if line.endswith('}'):
                self.finish_reaction()
                self.state = 0

        elif self.state == 2: # Continuation line
            #print("CONT", self.key, self.reaction.get(self.key, 'XXX'))
            if line.endswith('&'):
                self.reaction[self.key] = self.join(
                    self.reaction[self.key], line[:-1].strip())
                self.state = 2
            elif line.endswith('}'):
                if line != '}':
                    self.reaction[self.key] = self.join(
                        self.reaction[self.key], line[:-1].strip())
                self.finish_reaction()
                self.state = 0
            else:
                self.reaction[self.key] = self.join(
                    self.reaction[self.key], line)
                self.state = 1

    def join(self, l1, l2):
        #print("JOIN", l1, l2)
        # Fix adjoining quotes
        if l1[-1]==l2[0] and l2[0] in ('"',"'"):
            r = l1[:-1] + l2[1:]
        else:
            r =l1+l2
        #print("JOIN RETURNS", r)
        return r

    def finish_reaction(self):
        #print("FINISH", self.reaction_id)
        rxn = self.reactions[self.reaction_id]
        chem_eq = rxn.get('chem_eq','NONE')
        # Remove quotes
        if chem_eq[0]==chem_eq[-1] and chem_eq[0] in ('"',"'"):
            chem_eq = chem_eq[1:-1]
        rxn['reactants'], rxn['products'] = self.parse_chem_eq(chem_eq)
        if chem_eq.upper() == "NONE":
            rxn['chem_eq'] = "NONE"
        else:
            # Reformat chemical equation to canonicalize format
            fmt = {}
            for side in 'reactants', 'products':
                fmt[side] = ' + '.join(species if coeff==1.0 else '%g*%s' % (coeff, species)
                                       for (species, coeff) in rxn[side] if coeff)
            chem_eq = "%s --> %s" % (fmt['reactants'], fmt['products'])
            rxn['chem_eq'] =chem_eq

        for k, v in list(rxn.items()):
            if '(' in k:
                base_k = k.split('(')[0].strip()
                if base_k not in rxn:
                    rxn[base_k] = {}
                idx = k.split('(')[1].split(')')[0].strip()
                idx = int(idx)
                rxn[base_k][idx] = v
                del rxn[k]
        #print(self.reaction_id)
        #pprint(rxn)
        #print()
        self.reaction_id = None
        self.reaction = {}

    def is_alnum(self, s):
        return all((c.isalnum() or c=='_') for c in s)


    def valid_reaction_id(self, s):
        return len(s) <= 32 and self.is_alnum(s)


    def err(self, expected, got):
        raise ParseError('expected %s, got %s state=%s' % (expected, got, self.state))


    def parse_chem_eq(self, chem_eq):
        """parse chemical equation, return two lists (Reactants, Products)
        each list is a list of tuples (species, coefficient)"""
        # Remove quotes
        if chem_eq[0]==chem_eq[-1] and chem_eq[0] in ('"',"'"):
            chem_eq = chem_eq[1:-1]
        if chem_eq is None or chem_eq.upper()=='NONE':
            return [], []
        if '-->' in chem_eq:
            lhs, rhs = chem_eq.split('-->', 1)
        elif '==' in chem_eq:
            lhs, rhs = chem_eq.split('==', 1)
        else:
            raise ParseError('expected --> or == in chem_eq, got %s' % chem_eq)

        def parse_side(side):
            ret = []
            for field in side.split('+'):
                field = field.strip()
                if not field:
                    continue # ?
                if '*' in field:
                    coeff, species = field.split('*', 1)
                    ret.append([species.strip(), float(coeff)])
                elif ' ' in field:
                    coeff, species = field.split(' ', 1)
                    ret.append([species.strip(), float(coeff)])
                else:
                    coeff = ''
                    while field and field[0].isdigit() or field[0]=='.':
                        coeff += field[0]
                        field = field[1:]
                    if coeff:
                        coeff = float(coeff)
                    else:
                        coeff = 1.0
                    ret.append([field.strip(), coeff])
            return ret

        return parse_side(lhs), parse_side(rhs)


def main():
    from mfixgui.tools.reaction_parser_test_data import data
    rp = ReactionParser()
    rp.parse(data)
    return rp

if (__name__ == '__main__'):
    main()
