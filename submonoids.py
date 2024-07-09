from sys import argv
from itertools import combinations, chain, product

# from https://docs.python.org/3/library/itertools.html#itertools-recipes

def powerset(iterable):
    """Yield the powerset of an iterable."""
    s = list(iterable)
    for r in range(len(s) + 1):
        yield from map(set, combinations(s, r))

class IdempotentCommutativeMonoid:
    # we assume our monoids are commutative and idempotent, as is true of join lattices,
    # which allows faster checks of submonoidality
    def __init__(self, elements, mult, id):
        self.elements = elements
        self.mult = mult
        self.id = id

    def is_submonoid(self, subset) -> bool:
        """Check if a subset is a submonoid."""
        if not subset.issubset(self.elements):
            raise ValueError("given subset is not a subset of the monoid")

        if self.id not in subset:
            return False

        for (a, b) in combinations(subset, 2):
            if self.mult(a, b) not in subset:
                return False

        return True

    def enumerate_submonoids(self):
        """Enumerate a generator of submonoids."""
        for subset in powerset(self.elements):
            if self.is_submonoid(subset):
                yield subset

    def count_submonoids(self):
        """Count the number of submonoids."""
        return sum(1 for _ in self.enumerate_submonoids())

    def __mul__(self, other):
        """Compute the product monoid."""
        if not isinstance(other, IdempotentCommutativeMonoid):
            raise ValueError("Can only take the product of two monoids")

        elements = set(product(self.elements, other.elements))
        mult = lambda a, b: (self.mult(a[0], b[0]), other.mult(a[1], b[1]))
        id = (self.id, other.id)

        return IdempotentCommutativeMonoid(elements, mult, id)

def chain(n):
    """Yield a chain of n elements, viewed as a monoid under the join."""
    elements = set(range(n))
    mult = max
    id = 0
    return IdempotentCommutativeMonoid(elements, mult, id)

def lattice(dims):
    """Yield a lattice of the appropriate dimensions---note that an m x n lattice has
    (m + 1) x (n + 1) lattice points."""
    ret = None
    for dim in dims:
        if ret is None:
            ret = chain(dim + 1)
        else:
            ret *= chain(dim + 1)

    return ret

if __name__ == "__main__":
    print(lattice(map(int, argv[1:])).count_submonoids())
