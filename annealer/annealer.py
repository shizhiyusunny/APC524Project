class Annealer:
    def build(inputs):
        from . import constantStep
        if inputs['type'] == 'Constant Step':
            return constantStep.ConstantStep(inputs)
    def stepper(self):
        pass
    def store_all(processes, f, t=1.):
        for p in processes:
            p.store(f, t)
