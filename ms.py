import bisect
import heapq

import numpy as np

current_time = 0


class AncestralRecombinationGraph:
    global current_time

    def __init__(self, n, rho):
        """Initialize n sample lineages, which each span a recombination
        distance of rho."""
        global current_time
        current_time = 0
        self.lineages = Heap(
            Lineage(AncestralMaterial([{sample: [0, rho]}]))
            for sample in range(n))

    def simulate(self):
        global current_time
        current_time = 0

        n = len(self.lineages)
        while n > 1:
            coalescence_time = (current_time +
                np.random.exponential(scale=2.0 / n / (n - 1)))
            if coalescence_time < self.lineages[0].recombination_time:
                current_time = coalescence_time
                self.coalesce()
                event = 'coalescence'
            else:
                current_time = self.lineages[0].recombination_time
                self.recombine()
                event = 'recombination'
            n = len(self.lineages)
            print current_time, '\t', n, '\t', event

    def coalesce(self):
        index_a, index_b = np.random.choice(len(self.lineages), 2,
            replace=False)
        lineage_a = self.lineages[index_a]
        lineage_b = self.lineages[index_b]
        coalesced_lineage = Lineage.coalesce(lineage_a, lineage_b)

        del self.lineages[max(index_a, index_b)]
        del self.lineages[min(index_a, index_b)]
        self.lineages.heapify()
        self.lineages.insert(coalesced_lineage)

    def recombine(self):
        recombined_lineage = self.lineages.pop()
        lineage_a, lineage_b = recombined_lineage.recombine()
        self.lineages.insert(lineage_a)
        self.lineages.insert(lineage_b)


class Lineage:
    def __init__(self, ancestral_material):
        self.ancestral_material = ancestral_material
        self.recombination_time = (current_time +
            np.random.exponential(scale=1.0 / self.ancestral_material.span()))

    @staticmethod
    def coalesce(lineage_a, lineage_b):
        return Lineage(lineage_a.ancestral_material +
                       lineage_b.ancestral_material)

    def recombine(self):
        recombination_point = np.random.uniform(
            low=self.ancestral_material.inf,
            high=self.ancestral_material.sup)
        l_am, r_am = self.ancestral_material.split(recombination_point)

        return Lineage(l_am), Lineage(r_am)


class AncestralMaterial:
    def __init__(self, endpoints_maps, inf=None, sup=None):
        self.endpoints_maps = endpoints_maps

        self.inf = inf
        self.sup = sup
        if inf is None or sup is None:
            self.endpoints_maps = [self.endpoints_map()]
            self.inf = min(endpoints[0] for endpoints
                           in self.endpoints_maps[0].values())
            self.sup = max(endpoints[-1] for endpoints
                           in self.endpoints_maps[0].values())

    def __add__(self, other):
        return AncestralMaterial(self.endpoints_maps + other.endpoints_maps,
                                 inf=min(self.inf, other.inf),
                                 sup=max(self.sup, other.sup))

    def split(self, x):
        l_endpoints_map = {}
        r_endpoints_map = {}
        for sample, endpoints in self.endpoints_map().items():
            split_index = bisect.bisect(endpoints, x)

            l_endpoints = endpoints[:split_index]
            r_endpoints = endpoints[split_index:]
            if split_index % 2 == 1:
                l_endpoints = l_endpoints + [x]
                r_endpoints = [x] + r_endpoints

            if l_endpoints:
                l_endpoints_map[sample] = l_endpoints
            if r_endpoints:
                r_endpoints_map[sample] = r_endpoints

        return (AncestralMaterial([l_endpoints_map]),
                AncestralMaterial([r_endpoints_map]))

    def endpoints_map(self):
        def merge(endpoints_maps):
            if len(endpoints_maps) == 1:
                return endpoints_maps[0]
            else:
                midpoint = len(endpoints_maps) / 2
                return merge_2(merge(endpoints_maps[:midpoint]),
                               merge(endpoints_maps[midpoint:]))

        def merge_2(endpoints_map_a, endpoints_map_b):
            result = {}
            result.update(endpoints_map_a)
            result.update(endpoints_map_b)

            for sample in (set(endpoints_map_a.keys()) &
                          set(endpoints_map_b.keys())):
                result[sample] = merge_endpoints(endpoints_map_a[sample],
                                                 endpoints_map_b[sample])
            return result

        def merge_endpoints(endpoints_a, endpoints_b):
            """Merge two sorted lists into one sorted list, and discard elements common
            to both lists."""
            a = 0
            b = 0
            endpoints = []
            while (a < len(endpoints_a) and
                   b < len(endpoints_b)):
                if endpoints_a[a] < endpoints_b[b]:
                    endpoints.append(endpoints_a[a])
                    a += 1
                elif endpoints_a[a] > endpoints_b[b]:
                    endpoints.append(endpoints_b[b])
                    b += 1
                else:  # then the two elements are equal
                    a += 1
                    b += 1

            endpoints += endpoints_a[a:]
            endpoints += endpoints_b[b:]

            return endpoints

        return merge(self.endpoints_maps)

    def span(self):
        return self.sup - self.inf


class Heap:
    def __init__(self, lineages):
        self.items = [(lineage.recombination_time, lineage)
                      for lineage in lineages]
        self.heapify()

    def __len__(self):
        return len(self.items)

    def __delitem__(self, index):
        del self.items[index]

    def __getitem__(self, index):
        return self.items[index][1]

    def heapify(self):
        heapq.heapify(self.items)

    def pop(self):
        return heapq.heappop(self.items)[1]

    def insert(self, lineage):
        heapq.heappush(self.items, (lineage.recombination_time, lineage))


if __name__ == '__main__':
    arg = AncestralRecombinationGraph(1000, 10)
    arg.simulate()
