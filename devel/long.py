memo = {}

all_nodes = range(9)
from_nodes = ((), (0, ), (1, ), (2, 4), (0, ), (3, ), (5, 8), (3, ), (4, 7))


def get_longest(to_node):

    if to_node in memo:

        return memo[to_node]


    best = 0

    for from_node in from_nodes[to_node]:

        best = max(best, get_longest(from_node) + 1)

    memo[to_node] = best


    return best


length, node = max([(get_longest(to_node), to_node) for to_node in all_nodes])

print 'Length of longest path: length = %d, ending at %s' % (length, node)
