from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):
	def __init__(self, coordinates):
		try:
			if not coordinates:
				raise ValueError
			self.coordinates = tuple([Decimal(x) for x in coordinates])
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
		return Vector(new_coordinates)

	def minus(self,v):
		new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
		return Vector(new_coordinates)

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
			magnitude = self.magnitude()
			return self.times_scalar(Decimal('1.0')/Decimal(magnitude))
		except ZeroDivisionError:
			raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

	def dot_product(self, v):
		return sum([x*y for x,y in zip(self.coordinates, v.coordinates)])

	def angle(self, v,degrees=False):
		try:
			from math import pi,acos
			u = self.normalize()
			w = v.normalize()
			radians = acos(round(u.dot_product(w),10))

			if(degrees):
				degrees_per_radians = 180. / pi
				return radians * degrees_per_radians
			else:
				return radians

		except Exception as e:
			if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
				raise Exception('Cannot compute an angle with the zero vector')

	def isZeroVector(self):
		tolerance = 1e-10
		return (self.magnitude() < tolerance)

	def isParallel(self, v):
		if(self.isZeroVector() or v.isZeroVector()):
			return True
		scalar = self.coordinates[0] / v.coordinates[0]
		i = 1
		for i in range(len(self.coordinates)):
			if((round(scalar,10) != round(self.coordinates[i]/v.coordinates[i],13))):
				return False
		return True

	def isOrthogonal(self,v):
		tolerance = 1e-10
		return abs(self.dot_product(v)) < tolerance

	def component_parallelTo(self,basis):
		try:
			b = basis.normalize();
			dot = self.dot_product(b)
			return b.times_scalar(dot)
		except Exception as e:
			if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
				raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
			else:
				raise e

	def component_orthogonalTo(self,basis):
		try:
			return self.minus(self.projectionOn(basis))
		except Exception as e:
			if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
				raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
			else:
				raise e

	def crossProduct(self,v):
		try:
			from math import pi
			new_coordinates = []
			if(self.isZeroVector() or v.isZeroVector()):
				new_coordinates.append(0)
				new_coordinates.append(0)
				new_coordinates.append(0)		
			else:		
				new_coordinates.append((self.coordinates[1]*v.coordinates[2]) - (v.coordinates[1]*self.coordinates[2]))
				new_coordinates.append(-((self.coordinates[0]*v.coordinates[2]) - (v.coordinates[0]*self.coordinates[2])))
				new_coordinates.append((self.coordinates[0]*v.coordinates[1]) - (v.coordinates[0]*self.coordinates[1]))
			return Vector(new_coordinates)
		except ValueError as e:
			msg = str(e)
			if(msg == 'need more than 2 values to unpack'):
				self_embedded_in_R3 = Vector(self.coordinates + ('0',))
				v_embedded_in_R3 = Vector(v.coordinates + ('0',))
				return self_embedded_in_R3.crossProduct(v_embedded_in_R3)
			elif(msg == 'too many values to unpack' or msg == 'need more than 1 value to unpack'):
				raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
			else:
				raise e

	def parallelogramArea(self,v):
		vector = self.crossProduct(v)
		return vector.magnitude()

	def triangleArea(self,v):
		vector = self.crossProduct(v)
		return 0.5 * vector.magnitude()

v = Vector([1.5,9.547,3.691])
w = Vector([-6.007,0.124,5.772])
print(v.triangleArea(w))