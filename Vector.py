from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
	def __init__(self, coordinates):
		try:
			if not coordinates:
				raise ValueError
			self.coordinates = tuple(Decimal(x) for x in coordinates)
			self.dimension = len(coordinates)

		except ValueError:
			raise ValueError('The coordinates must be nonempty')

		except TypeError:
			raise TypeError('The coordinates must be an iterable')


	def __str__(self):
		return 'Vector: {}'.format(self.coordinates)


	def __eq__(self, v):
		return self.coordinates == v.coordinates

	def plus(self,v):
		new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordenates)

	def minus(self,v):
		new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordenates)

	def times_scalar(self,scalar):
		new_coordinates = [x * Decimal(scalar) for x in self.coordinates]
		return Vector(new_coordinates)

	def magnitude(self):
		acc = 0
		for i in range(len(self.coordinates)):    		
			acc += self.coordinates[i]**2   		
		from math import sqrt
		return sqrt(acc)

	def normalize(self):
		try:
			return self.times_scalar(Decimal('1.0')/magnitude)
		except ZeroDivisionError:
			raise Exception('Cannot normalize zero vector')

	def dot_product(self, v):
		acc = 0
		for i in range(len(self.coordinates)):
			acc = acc + self.coordinates[i] * v.coordinates[i]
		return acc

	def angle(self, v,degrees=True):
		try:
			from math import acos, pi
			if degrees:
				degrees = 180 / pi
				return acos(self.dot_product(v)/(self.magnitude() * v.magnitude())) * degrees
			else:
				return acos(self.dot_product(v)/(self.magnitude() * v.magnitude()))

		except ZeroDivisionError:
			raise Exception('Cannot compute an angle with a zero vector')

v = Vector([7.35,0.221,5.188])
w = Vector([2.751,8.259,3.985])
print (v.angle(w))
