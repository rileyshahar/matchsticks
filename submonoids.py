from sys import argv
from itertools import combinations, chain, product
from functools import lru_cache

# from https://docs.python.org/3/library/itertools.html#itertools-recipes

def powerset(iterable):
    """Yield the powerset of an iterable."""
    s = list(iterable)
    for r in range(len(s) + 1):
        yield from map(frozenset, combinations(s, r))

class IdempotentCommutativeMonoid:
    # we assume our monoids are commutative and idempotent, as is true of join lattices,
    # which allows faster checks of submonoidality
    def __init__(self, elements, mult, id):
        self.elements = elements
        self.mult = mult
        self.id = id

    @lru_cache
    def is_submagma(self, subset) -> bool:
        """Check if a subset is a submagma."""
        if not subset.issubset(self.elements):
            raise ValueError("given subset is not a subset of the monoid")

        for (a, b) in combinations(subset, 2):
            if self.mult(a, b) not in subset:
                return False

        return True

    @lru_cache
    def is_submonoid(self, subset) -> bool:
        """Check if a subset is a submonoid."""
        if self.id not in subset:
            return False

        if not self.is_submagma(subset):
            return False

        return True

    @lru_cache
    def is_ideal(self, subset) -> bool:
        """Check if a subset is an ideal."""
        if not subset.issubset(self.elements):
            raise ValueError("given subset is not a subset of the monoid")

        for a in self.elements:
            for b in subset:
                if self.mult(a, b) not in subset:
                    return False

        return True

    def enumerate_submagmas(self):
        """Enumerate a generator of submagmas."""
        for subset in powerset(self.elements):
            if self.is_submagma(subset):
                yield subset

    def enumerate_submonoids(self):
        """Enumerate a generator of submonoids."""
        for subset in powerset(self.elements):
            if self.is_submonoid(subset):
                yield subset

    def enumerate_ideals(self):
        """Enumerate a generator of ideals."""
        for subset in powerset(self.elements):
            if self.is_ideal(subset):
                yield subset

    def count_submagmas(self):
        """Count the number of submagmas."""
        return sum(1 for _ in self.enumerate_submagmas())

    def count_submonoids(self):
        """Count the number of submonoids."""
        return sum(1 for _ in self.enumerate_submonoids())

    def grade_sequence(self) -> list[int]:
        """Return the count of submonoids of each grade."""
        counts = [0 for _ in range(len(self.elements) + 1)]
        for submagma in self.enumerate_submonoids():
            counts[len(submagma)] += 1

        return counts

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

def partitioned_sum(m: int, p: IdempotentCommutativeMonoid):
    ret = 0
    for a in lattice.enumerate_submagmas():
        ret += c(m, p, IdempotentCommutativeMonoid(a, p.mult, None))

@lru_cache
def transition(a: IdempotentCommutativeMonoid, b: frozenset) -> int:
    """B is a submagma of A."""

    ret = 0
    for u in a.enumerate_ideals():
        if u.union(b) == a.elements:
            ret += 1

    return ret

def make_matrix(lattice: IdempotentCommutativeMonoid):
    return [
        [transition(
            IdempotentCommutativeMonoid(a, lattice.mult, None),
            b
        ) for a in lattice.enumerate_submonoids()]
    for b in lattice.enumerate_submonoids()]

def make_matrix_from_dims(dims: list):
    return make_matrix(lattice(dims))

def sub_nx1(k: int, S: frozenset, l: int):
    counter = 0
    elements = set()
    while counter < l:
        elements.add((counter, 0))
        counter += 1

    while counter < k + l:
        elements.add((counter, 1))
        if counter - l in S:
            elements.add((counter, 0))
        counter += 1

    def pair_max(a, b):
        return (max(a[0], b[0]), max(a[1], b[1]))

    return IdempotentCommutativeMonoid(frozenset(elements), pair_max, None)

def print_nx1(elements):
    if not elements:
        return "∅"
    l = max(a for a, _ in elements)
    s = [" " for _ in range (l * 2 + 1)]
    s2 = [" " for _ in range (l * 2 + 1)]
    for a, b in elements:
        if b == 1:
            s[a * 2] = "·"
        else:
            s2[a * 2] = "·"

    return "".join(s) + "\n" + "".join(s2)

def iso_class(elements):
    if not elements:
        return (0, frozenset(), 0)

    k_set = {a for a, b in elements if b == 1}
    s_set = {a for a, b in elements if b == 0 and a in k_set}
    l_set = {a for a, b in elements if b == 0 and a not in k_set}
    l = len(l_set)
    k = len(k_set)
    s = frozenset(a - min(k_set) for a in s_set) if k_set else frozenset()

    return (k, s, l)



def iso_transition(a: IdempotentCommutativeMonoid, k, S, l) -> int:
    """B is a submagma of A."""
    ret = 0
    for b in a.enumerate_submagmas():
        if (k, S, l) == iso_class(b):
            for u in a.enumerate_ideals():
                if u.union(b) == a.elements:
                    ret += 1
    return ret

def enumerate_subisos(k, S, l):
    for k2 in range(k + 1):
        for S2 in powerset(range(k2)):
            for l2 in range(k + len(S) + 1):
                yield (k2, S2, l2)

if __name__ == "__main__":
    a = sub_nx1(4, {0, 3}, 2)
    print(print_nx1(a.elements))
    for iso in enumerate_subisos(*iso_class(a.elements)):
        if (i := iso_transition(a, *iso)) > 0:
            print("------")
            print(iso[0], set(iso[1]), iso[2])
            print(i)
