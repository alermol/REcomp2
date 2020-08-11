class QuickUnion:

    def __init__(self, N):
        self._id = list(range(N))
        self._count = N
        self._rank = [0] * N

    def find(self, p):
        id = self._id
        while p != id[p]:
            p = id[p] = id[id[p]]   # Path compression using halving.
        return p

    def union(self, p, q):
        id = self._id
        rank = self._rank
        i = self.find(p)
        j = self.find(q)
        if i == j:
            return
        self._count -= 1
        if rank[i] < rank[j]:
            id[i] = j
        elif rank[i] > rank[j]:
            id[j] = i
        else:
            id[j] = i
            rank[i] += 1

    def __str__(self):
        return " ".join([str(x) for x in self._id])
