class Annealer:
    def build(inputs):
        from . import constantStep
        if inputs['type'] == 'Constant Step':
            return constantStep.ConstantStep(inputs)
    def stepper():
        pass
