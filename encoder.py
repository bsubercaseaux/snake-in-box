from eznf import modeler
import itertools
import argparse

def encode(n, k):
    enc = modeler.Modeler()

    vertices = list(itertools.product([0, 1], repeat=n))
    graph = {v :[] for v in vertices}
    for v in vertices:
        for i in range(n):
            neighbor = list(v)
            neighbor[i] = 1 - neighbor[i]
            graph[v].append(tuple(neighbor))


    # for u in vertices:
    #     for v in graph[u]:
    #         if u < v:
    #             enc.add_var(f"e_{u, v}")

    for v in vertices:
        for i in range(k):
            enc.add_var(f"x_{v, i}")
            enc.add_var(f"xle_{v, i}")

    for v in vertices:
        for i in range(k):
            enc.add_clause([f"-x_{v, i}", f"xle_{v, i}"])
            if i < k-1:
                enc.add_clause([f"-x_{v, i}", f"-xle_{v, i+1}"])
            if i > 0:
                enc.add_clause([f"-xle_{v, i}", f"xle_{v, i-1}"])    
        

    for v in vertices:
        enc.at_most_one([f"x_{v, i}" for i in range(k)])


    for i in range(k):
        enc.exactly_one([f"x_{v, i}" for v in vertices])

    

    for v in vertices:
        for i in range(1, k):
            enc.add_clause([f"-x_{v, i}"] + [f"x_{u, i-1}" for u in graph[v]])
    
            for u in graph[v]:
                enc.add_clause([f"-x_{v, i}", f"xle_{u, i-1}"])
        # for i in range(2, k):
        #     for u in graph[v]:
        #         for j in range(i-1):
        #             enc.add_clause([f"-x_{v, i}", f"-x_{u, j}"])

    # symmetry breaking
    for i in range(n):
        for v in vertices:
            for j in range(i, n):
                if v[j] == 1:
                    enc.add_clause([f"-x_{v, i}"])



    return enc


def decode(model, n, k):
    vertices = list(itertools.product([0, 1], repeat=n))
  
    sequence = []
    for i in range(k):
        for v in vertices:
            if model[f"x_{v, i}"]:
                sequence.append(v)

    print(f"Decoded sequence: {sequence}")
    assert len(sequence) == k, f"Decoded sequence length {len(sequence)} does not match {k}"
    return sequence

argparser = argparse.ArgumentParser(description="Snake Box Encoder")
argparser.add_argument("-n", type=int, help="Number of vertices")
argparser.add_argument("-k", type=int, help="Number of positions")
argparser.add_argument("--solve", action="store_true", help="Solve the encoding")
args = argparser.parse_args()
n = args.n
k = args.k
encoding = encode(n, k)
if args.solve:
    encoding.solve_and_decode(lambda model: decode(model, n, k))
else:
    encoding.serialize(f"snake-box/s_{n}_{k}.cnf")
