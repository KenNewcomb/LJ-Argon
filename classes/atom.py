# An object that represents a single "atom" to be studied.
# Each atom has Cartesian coordinates x, y, and z, as well as their respective
# components of velocity (vx, vy, and vz).

class Atom:
	# 6 variables: 3 cartestian coordinates, 3 components of velocity.
	x = 0
	y = 0
	z = 0
	vx = 0
	vy = 0
	vz = 0

	def __init__(self):
		"""Constructs an atom"""
	
	def randcoord(self):
		pass
	# def the_boltzmannator(self):
