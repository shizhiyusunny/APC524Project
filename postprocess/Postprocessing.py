class Postprocessing:
    def build(inputs):
        processes = []
        from . import storePVD,storeNPY,storePLT
        for i in inputs:
            if i['type'] == 'pvd':
                processes.append(storePVD.StorePVD(i))
            elif i['type']=="plot":
                processes.append(storePLT.StorePLT(i))
        return processes
    def __init__(self, inputs):
        self.iterative = inputs['iterative']
    def store(self, f, t):
        pass
