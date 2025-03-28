
import json, re
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

    form_name_map = dict(map(str.split,filter(lambda s: s.strip(),"""
        Fucp	Fuc
        Fucx	Fuc
        Galf	Gal
        Galp	Gal
        GalpNAc	GalNAc
        GalxNAc	GalNAc
        GlcAp	GalA
        Glcp	Glc
        GlcpNAc	GlcNAc
        Glcx	Glc
        GlcxNAc	GlcNAc
        KDNp	KDN
        Manp	Man
        ManpNAc	ManNAc
        Manx	Man
        NeupNAc NeuAc
        NeupNGc NeuGc
        Xylf	Xyl
        Xylp	Xyl
    """.splitlines())))

    def map_form_names(self,row):
        for fn in ("form_name","parent_form_name"):
            row[fn] = self.form_name_map.get(row[fn],row[fn])
        return row
    
    def group(self,enzyme):
        return self.groupby%enzyme
        
    def __init__(self, datafile=None, groupby=None, genesets=None):
        if not datafile:
            self.datafile = self.sandbox_datafile
        self.datafile = datafile
        self.enzyme2group = dict()
        self.group2enzymes = defaultdict(set)
        self.groupby = self.convert_group(groupby)
        self.genesets_generators = [ self.genesets_generator_map[gs.strip()] for gs in genesets.split(',') ]
        self.build()

    def build(self):

        h = open(self.sandbox_datafile)
        data = json.loads(h.read())
        h.close()
        
        byid = defaultdict(list)
        for row in data["data"]:
            if row['species'] != 'Homo sapiens' or row['gene_name'] == None:
                continue
            row = self.map_form_names(row)
            byid[row['residue_id']].append(dict(row.items()))
            group = self.group(row)
            self.enzyme2group[row['gene_name']] = group
            self.group2enzymes[group].add(row['gene_name'])
        self.adjgroups = set()
        for rid in byid:
            for row in byid[rid]:
                if row['parent_id'] == "no_id":
                    continue
                for prow in byid.get(row['parent_id'],[]):
                    group = self.group(row)
                    pgroup = self.group(prow)
                    self.adjgroups.add((group,pgroup))

    def all_enzymes(self):
        return list(self.enzyme2group)

    def all_groups(self):
        return list(self.group2enzymes)

    def enzymes_bygroup(self,group):
        return list(self.group2enzymes.get(group,[]))

    def common_prefix_genesetname(self,geneset):
        if len(geneset) <= 1:
            return ",".join(geneset)
        genelist = sorted(geneset)
        valid_pflen = []
        for gene in genelist:
            m = re.search(r"\d+$",gene)
            if m:
                 valid_pflen.append(len(gene)-len(m.group(0))+1)
            else:
                 valid_pflen.append(len(gene))
        pflen = 0
        for l in range(1,min(valid_pflen)):
            good = True
            for i in range(1,len(genelist)):
                if genelist[i][:l] != genelist[0][:l]:
                    good = False
                    break
            if good:
                pflen = l
            else:
                break
        if pflen == 0:
            return ",".join(genelist)
        remainders = [ gene[l:] for gene in genelist ]
        remasint = sorted(map(int,[ r for r in remainders if re.search(r'^\d+$',r) ]))
        remasstr = sorted([ r for r in remainders if not re.search(r'^\d+$',r) ])
        if len(remasint) > 1:
            remlst = []
            while len(remasint) > 0:
                ed = 0
                for i in range(1,len(remasint)):
                    if remasint[ed] != (remasint[i]-1):
                        break
                    ed = i
                if ed == 0:
                    remlst.append(str(remasint[0]))
                else:
                    remlst.append("%s-%s"%(remasint[0],remasint[ed]))
                remasint = remasint[(ed+1):]
            remstr = ",".join(remlst) + "," + ",".join(remasstr)
            remstr = remstr.strip(",")
        else:
            remstr = ",".join(map(str,remasint)) + "," + ",".join(sorted(remasstr))
            remstr = remstr.strip(",")
        genesetname = genelist[0][:l] + "{" + remstr + "}"
        return genesetname

    def genesets(self):
        for gsgen in self.genesets_generators:
            for gsname,gs in gsgen(self):
                yield gsname,gs

    def singletons(self):
        for grp in self.all_groups():
            for gn in self.enzymes_bygroup(grp)
                yield grp,set([gn])

    def pairs(self):
        for gn1 in self.all_enzymes():
            for gn2 in self.all_enzymes():
                if gn1 < gn2:
                    enzs = set([gn1,gn2])
                    yield self.common_prefix_genesetname(enzs),enzs

    def groupsets(self):
        for grp in self.all_groups():
            enzs = self.enzymes_bygroup(grp)
            yield grp,enzs

    def nongrouppairs(self):
        for grp1 in self.all_groups():
            for grp2 in self.all_groups():
                if grp1 >= grp2:
                    continue
                for gn1 in self.enzymes_bygroup(grp1):
                    for gn2 in self.enzymes_bygroup(grp2):
                        enzs = set([gn1,gn2])
                        yield ",".join([grp1,grp2]),enzs

    def grouppairs(self):
        for grp in self.all_groups():
            for gn1 in self.enzymes_bygroup(grp):
                for gn2 in self.enzymes_bygroup(grp):
                    if gn1 < gn2:
                        enzs = set([gn1,gn2])
                        yield grp,enzs

    def adjgrouppairs(self):
        for grp1,grp2 in self.adjgroups:
            for gn1 in self.enzymes_bygroup(grp1):
                for gn2 in self.enzymes_bygroup(grp2):
                    enzs = set([gn1,gn2])
                    yield "<-".join([grp1,grp2]),enzs

    genesets_generator_map = dict(
        SINGLETONS = singletons,
        NONGROUPPAIRS = nongrouppairs,
        ADJGROUPPAIRS = adjgrouppairs,
        GROUPPAIRS = grouppairs,
        PAIRS = pairs,
        BYGROUP = groupsets,
    )


if __name__ == "__main__":
    ge = GlycoEnzymes("data/sandbox.json","GROUP_MONO_ANOMER_SITE_PARENT","BYGROUP")
    print(ge.all_groups())
    print(ge.adjgroups)

