from pysat.pb import *
from pysat.formula import CNF, IDPool
from pysat.solvers import Cadical195, Glucose42
from pysat.card import CardEnc
import itertools
import argparse
from lex import lex_smaller_eq, checkLexMin
# from checker import check_hamiltonian_path



def permute(v, i, j):
        v_list = list(v)
        v_list[i], v_list[j] = v_list[j], v_list[i]
        return tuple(v_list)
    
def flip_i(v, i):
    v_list = list(v)
    v_list[i] = 1 - v_list[i]
    return tuple(v_list)


def at_most_one(variables, vpool):
    if len(variables) <= 4:
        return [[-comb[0], -comb[1]] for comb in itertools.combinations(variables, 2)]
    else:
        new_var = vpool.id()
        return at_most_one(variables[:3] + [-new_var], vpool) + at_most_one([new_var] + variables[3:], vpool)
    

def base_induced_path(G, s, k):
    """
    Base encoding for Induced path in a graph G from node s.

    :param G: The graph represented as an adjacency list.
    :param s: The starting node of the path.
    :param k: The length of the path.
    """
    cnf = CNF()
    n = len(list(G.keys())[0])

    vertices = list(sorted(G.keys()))
    
    # Create edge variables mapping
    edges = []
    edge_vars = {}
    vpool = IDPool()
    for u in vertices:
        for v in G[u]:
            edges.append((u, v))

    ev = lambda e: vpool.id(f"e_{e[0], e[1]}")

    p = lambda v: vpool.id(f"p_{v}")

    p_pos = lambda v, i: vpool.id(f"p_{v}_{i}")

    last = lambda v: vpool.id(f"last_{v}")

    for e in edges:
        cnf.append([-ev(e), +p(e[0])])  
        cnf.append([-ev(e), +p(e[1])])

    cnf.extend(CardEnc.equals(lits=[last(v) for v in vertices], bound=1, vpool=vpool).clauses)  # at least one vertex is present
    # exact number of edges
    for u in vertices:
        cnf.extend(at_most_one([-p(u)] + [ev((u, v)) for v in G[u]], vpool)) # <= 1 outgoing edge
        cnf.extend(at_most_one([-p(u)] + [ev((v, u)) for v in G[u]], vpool)) # <= 1 incoming edge
        if u != s:
            cnf.append([-p(u), +last(u)] + [ev((u, v)) for v in G[u]]) # >= 1 outgoing edge
            cnf.append([-p(u)] + [ev((v, u)) for v in G[u]]) # >= 1 incoming edge
          
        if u == s:
            cnf.append([ev((s, v)) for v in G[s]]) # >= 1 outgoing edge
            cnf.extend([[-ev((v, s))] for v in G[s]])  # no incoming edges to s
        

        for v in G[u]:
            cnf.append([-last(u), -ev((u, v))])

    # no backtracks
    for u in vertices:
        for v in G[u]:
            cnf.append([-ev((u, v)), -ev((v, u))])  # if u -> v then not v -> u

    cnf.append([p_pos(s, 0)])
    for v in vertices:
        for u in G[v]:
            for i in range(n-1):
                cnf.append([-p_pos(v, i), -ev((v, u)), p_pos(u, i+1)])  # if v is at position i and v -> u then u is at position i+1

    for i in range(n):
        cnf.extend(CardEnc.equals(lits=[p_pos(v, i) for v in vertices], bound=1, vpool=vpool).clauses)  # 



    cnf.append([p(s)])  # starting vertex must be present
    cnf.extend(CardEnc.atleast(lits=[p(v) for v in vertices], bound=k+1, vpool=vpool, encoding=7).clauses)  

    edge_vars = [ev((u, v)) for u in vertices for v in G[u]]
    cnf.extend(CardEnc.atleast(lits=edge_vars, bound=k, vpool=vpool, encoding=7).clauses)

    for i in range(n):
        cnf.append([p(v) for v in vertices if v[i] == 1])  # at least one vertex must have bit i set to 1
    # no possible shortcut extensions
    # for u in vertices:
    #     for v in G[u]:
    #         for i in range(n):
    #             if u[i] != v[i]: 
    #                 continue
    #             up, vp = flip_i(u, i), flip_i(v, i)
    #             cls = [p(up), p(vp), -ev((u, v))]
    #             for x in G[vp]:
    #                 if x != v:
    #                     cls.append(p(x))
    #             cnf.append(cls)

    for u in vertices:
        for v in G[u]:
            if u < v:
                cnf.append([-p(v), -p(u), ev((u, v)), ev((v, u))])  # if p(v) and p(u) then e(u, v)

    # symmetry breaking
    for i in range(n):
        for v in vertices:
            for j in range(n-i):
                if v[j] == 1:
                    cnf.append([-p_pos(v, i)]) 
    # original_vertices = [v for v in vertices]
    # for i, j in itertools.combinations(range(n), 2):
    #     # print(f"Permuting dimensions {i} and {j} ", len(v))
    #     permuted_vertices = [permute(v, i, j) for v in original_vertices]
    #     lex_smaller_eq(cnf, vpool, [p(v) for v in original_vertices], [p(v) for v in permuted_vertices])
    # print(f"len(cnf) = {len(cnf.clauses)}")
    # for i in range(n):
    #     permuted_vertices = [flip_i(v, i) for v in original_vertices]
    #     lex_smaller_eq(cnf, vpool, [p(v) for v in original_vertices], [p(v) for v in permuted_vertices])

    return cnf, vpool

def incremental_induced_path(G, s, k, verbose=False):
    base_cnf, vpool = base_induced_path(G, s, k)
    vertices = list(sorted(G.keys()))
    solver = Cadical195(base_cnf)
    # solver = Glucose42(base_cnf)
    formula = base_cnf

    cycle_counter = 0
    round_counter = 0

    total_cycles = []
    found = False

    while True:
        if verbose:
            n_clauses = solver.nof_clauses()
            print(n_clauses, "clauses in the solver")
            if n_clauses > 100000:
                formula.to_file(f"induced_path_{n_clauses}.cnf")
            # if solver.nof_clauses() > 130000:
            #     solver.

        solver.solve()
        model = solver.get_model()
        if model is None:
            print("No induced path found.")
            break 
        
        else:
            # print(f"Model found with {len(model)} variables.")
            # print(f"Variables in the model: {model}")
            edges = []
            for u in vertices:
                for v in G[u]:
                    # print(f"Checking edge {u} -- {v}: id = {vpool.id(f'e_{u, v}')}")
                    if model[vpool.id(f"e_{u, v}") - 1] > 0:
                        # print(f"Edge {u} -- {v} is part of the Hamiltonian path.")
                        edges.append((u, v))
            # identify cycles now
            path = [s]
            current = s
            visited = {s}
            while True:
                found = False
                for v in G[current]:
                    if (current, v) in edges or (v, current) in edges:
                        if v not in visited:
                            path.append(v)
                            visited.add(v)
                            current = v
                            found = True
                            break
                if not found:
                    if verbose:
                        print("No further path found, terminating early.")
                        print(f"Vertex {current} has no unvisited neighbors.")
                    break
            # check length
            if len(path) - 1 >= k :
                if verbose:
                    print(f"Induced path found: {' -> '.join(map(str, path))}")
                    print(f"present edges = {edges}")
                    found = True
                    # print(f"Checking induced path validity...")
                    break
                    # return path
            else:
                print(f"Induced path found is not of length {k}, there must be a cycle")
                print(f"Path length: {len(path)-1}, expected: {k}")
                print(f"Path: {' -> '.join(map(str, path))}")
                round_counter += 1
                cycles = []
                for node in G:
                    if node not in visited:
                        cycle = [node]
                        current = node
                        while True:
                            found = False
                            for v in G[current]:
                                if (current, v) in edges:
                                    if v not in cycle:
                                        cycle.append(v)
                                        current = v
                                        found = True
                                        visited.add(v)
                                        break
                            if not found:
                                if len(cycle) == 1: break
                                cycles.append(cycle)
                                if verbose:
                                    print(f"Cycle found: {' -> '.join(map(str, cycle))}")
                                    print(f"So far cycle_counter = {cycle_counter}, and round_counter = {round_counter}")
                                solver.add_clause([-vpool.id(f"e_{cycle[i], cycle[i+1]}") for i in range(len(cycle)-1)] + [-vpool.id(f"e_{cycle[-1], cycle[0]}")])
                                formula.append([-vpool.id(f"e_{cycle[i], cycle[i+1]}") for i in range(len(cycle)-1)] + [-vpool.id(f"e_{cycle[-1], cycle[0]}")])
                                cycle_counter += 1
                                break
                total_cycles.extend(cycles)
    print(total_cycles)
    cycle_lengths = [len(c) for c in total_cycles]
    min_cycle_length = min(cycle_lengths) if cycle_lengths else 0
    max_cycle_length = max(cycle_lengths) if cycle_lengths else 0
    avg_cycle_length = sum(cycle_lengths) / len(cycle_lengths) if cycle_lengths else 0
    print(f"Total cycles found: {len(total_cycles)}")
    print(f"Min cycle length: {min_cycle_length}")
    print(f"Max cycle length: {max_cycle_length}")
    print(f"Avg cycle length: {avg_cycle_length:.2f}")
    # print(f"Cycle lengths: {cycle_lengths}")
    # print(f"Avg cycle length = {sum(len(c) for c in total_cycles)/len(total_cycles) if total_cycles else 0}")

    if found:
        return path
    else:
        return None


def main():
    argparser = argparse.ArgumentParser(description="Snake Box Encoder")
    argparser.add_argument("-n", type=int, help="Number of vertices", required=True)
    argparser.add_argument("-k", type=int, help="Number of positions", required=True)
    argparser.add_argument("--solve", action="store_true", help="Solve the encoding")
    args = argparser.parse_args()
    n = args.n
    k = args.k
    vertices = list(itertools.product([0, 1], repeat=n))
    graph = {v :[] for v in vertices}
    for v in vertices:
        for i in range(n):
            neighbor = list(v)
            neighbor[i] = 1 - neighbor[i]
            graph[v].append(tuple(neighbor))

    incremental_induced_path(graph, vertices[0], k, verbose=True)

if __name__ == "__main__":
    main()