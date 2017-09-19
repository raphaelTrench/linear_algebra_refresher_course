from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = system.dimension

        col = 0
        for row in range(num_equations):
            while col < num_variables:
                c = MyDecimal(system[row].normal_vector[col])
                if c.is_near_zero():
                    swap_succeeded = system.swapWithNonZeroBelow(row, col)
                    if not swap_succeeded:
                        col += 1
                        continue

                system.clearCoefficientsBelow(row, col)
                col += 1
                break

        return system

    def compute_rref(self):
        tf = self.compute_triangular_form()
        index = tf.indices_of_first_nonzero_terms_in_each_row()
        num_equations = len(tf)        

        for i in range(num_equations)[::-1]:
            j = index[i]
            if (j<0):
                continue
            tf.scaleRowToMakeCoefficientOne(i, j)
            tf.clearCoefficientsAbove(i, j)
        return tf

    def scaleRowToMakeCoefficientOne(self,row,col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)

    def clearCoefficientsAbove(self,row,col):
        for row_to_be_added_to in range(row)[::-1]:
            n = self[row_to_be_added_to].normal_vector
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, row_to_be_added_to)

    def swapWithNonZeroBelow(self,row,col):
        num_equations = len(self)

        for k in range(row + 1, num_equations):
            coefficient = MyDecimal(self[k].normal_vector[col])
            if not coefficient.is_near_zero():
                self.swap_rows(row, k)
                return True
        return False    

    def clearCoefficientsBelow(self,row,col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector[col])

        for row_to_be_added_to in range(row + 1, num_equations):
            n = self[row_to_be_added_to].normal_vector
            gamma = n[col]
            alpha = -gamma / beta
            self.add_multiple_times_row_to_row(alpha, row, row_to_be_added_to)

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        aux = self[row].normal_vector.times_scalar(coefficient)
        aux2 = self[row].constant_term * coefficient
        self[row] = Plane(normal_vector=aux, constant_term=aux2)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        row1 = self[row_to_add].normal_vector.times_scalar(coefficient)
        k1 = self[row_to_add].constant_term * coefficient
        row2 = self[row_to_be_added_to].normal_vector.plus(row1)
        k2 = self[row_to_be_added_to].constant_term
        k = k1 + k2
        self[row_to_be_added_to] = Plane(normal_vector=row2,constant_term=k)        

   

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='-1') and
        r[1] == p2):
    print ('test case 1 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print ('test case 2 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        r[3] == Plane()):
    print ('test case 3 failed')

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
    print ('test case 4 failed')
