import cogent.db.ensembl
class GeneRecord():
    
    def __init__(self, header, seq, count):
        self.header = header.rstrip()
        self.seq = seq.replace("\n","")
        self.count = count

    def __str__(self):
        return self.header+"| match count:"+str(self.count)+"\n"+self.seq

    def no_seq(self):
        return self.header+"\t"+str(self.count)