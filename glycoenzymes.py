
import json
from collections import defaultdict


class GlycoEnzymes():

    sandbox_datafile = "data/sandbox.json"
    
    GROUP_MONO_ANOMER_SITE_PARENT = "%(form_name)s-%(anomer)s-%(site)s-%(parent_form_name)s"
    GROUP_MONO_ANOMER_PARENT = "%(form_name)s-%(anomer)s-%(parent_form_name)s"
    GROUP_MONO_SITE_PARENT = "%(form_name)s-%(site)s-%(parent_form_name)s"
    GROUP_MONO_PARENT = "%(form_name)s-%(parent_form_name)s"
    GROUP_MONO_ANOMER_SITE = "%(form_name)s-%(anomer)s-%(site)s"
    GROUP_MONO_ANOMER = "%(form_name)s-%(anomer)s"
    GROUP_MONO_SITE = "%(form_name)s-%(site)s"
    GROUP_MONO = "%(form_name)s"

    @staticmethod
    def convert_group(s):
        return getattr(GlycoEnzymes,s)
    
    def group(self,enzyme):
        return self.groupby%enzyme
        
    def __init__(self, datafile=None, groupby=None):
        if not datafile:
            self.datafile = self.sandbox_datafile
        self.datafile = datafile
        self.enzyme2group = dict()
        self.group2enzymes = defaultdict(set)
        self.groupby = self.convert_group(groupby)

        self.build()

    def build(self):

        h = open(self.sandbox_datafile)
        data = json.loads(h.read())
        h.close()
        
        for row in data["data"]:
            if row['species'] != 'Homo sapiens' or row['gene_name'] == None:
                continue
            group = self.group(row)
            self.enzyme2group[row['gene_name']] = group
            self.group2enzymes[group].add(row['gene_name'])

    def all_enzymes(self):
        return list(self.enzyme2group)

    def all_groups(self):
        return list(self.group2enzymes)

    def enzymes_bygroup(self,group):
        return list(self.group2enzymes.get(group,[]))

    def singletons(self):
        for gn in self.all_enzymes():
            yield gn,set([gn])

    def pairs(self):
        for gn1 in self.all_enzymes():
            for gn2 in self.all_enzymes():
                if gn1 < gn2:
                    yield "%s,%s"%(gn1,gn2),set([gn1,gn2])

    def groupsets(self):
        for grp in self.all_groups():
            yield grp,self.enzymes_bygroup(grp)

