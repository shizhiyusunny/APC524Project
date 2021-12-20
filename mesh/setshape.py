from fenics import *
from mshr import *


class SetShape:

    def __init__(self, domain1, domain2):
        self.domain1 = domain1
        self.domain2 = domain2

    def addiction(self):
        domain = self.domain1 + self.domain2
        return domain

    def subtraction(self):
        domain = self.domain1 - self.domain2
        return domain

