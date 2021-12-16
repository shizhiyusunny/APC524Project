class Shape:
    def build(inputs):
        from . import rectangle, circle
        if 'Rectangle' in inputs:
            return rectangle.Rectangle(inputs['Rectangle'])
        elif 'Circle' in inputs:
            return circle.Circle(inputs['Circle'])
