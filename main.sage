import bisect;

p = 7;
F = FiniteField(p);

def add(a, b):
    return [a[i] + b[i] for i in range(len(a))];

def add_set(a, X):
    return [[a[i] + b[i] for i in range(len(a))] for b in X];

def sub(a, b):
    return [a[i] - b[i] for i in range(len(a))];

def dot(a, b):
    s = 0;
    for i in range(len(a)):
        s = s + a[i] * b[i];
    return s;
;

def quadrance(A, B):
    a = sub(B, A);
    return dot(a, a);
;

def angle(A, O, B):
    a = sub(A, O);
    b = sub(B, O);
    return dot(a, b); # / (a * b)

def unit_circle():
    A = [];
    for i in F:
        if (1 + i^2) != 0: # tangent
            c = (1 - i^2)/(1 + i^2); # cosine
            s = 2*i / (1 + i^2); # sine
            A.append([c, s]);
    A.append([-1 * (F.from_integer(1)), F.from_integer(0)]);
    return A;
;

def octant_sphere():
    H = []; # half of field
    for i in F:
        H.append(i);
    i = 0;
    while i < len(H):
        H = subtract_(H, [-H[i]]);
        i = i + 1;
    H.append(0);
    S = [];
    for i in H:
        for j in H:
            for k in H:
                if i^2 + j^2 + k^2 == 1:
                    S.append([i, j, k]);
    return S;
;

def unit_sphere(S):
    T = [];
    for i in S:
        T.append(i);
        T.append([-i[0], -i[1], -i[2]]);
        a0 = (i[0] == 0);
        b0 = (i[0] == 0);
        c0 = (i[0] == 0);
        ; if a0 else T.append([-i[0], i[1], i[2]]);
        ; if a1 else T.append([i[0], -i[1], i[2]]);
        ; if a2 else: T.append([i[0], i[1], -i[2]]);
        ; if a0 and a1 else: T.append([-i[0], -i[1], i[2]]);
        ; if a0 and a2 else: T.append([-i[0], i[1], -i[2]]);
        ; if a1 and a2 else: T.append([i[0], -i[1], -i[2]]);
    return T;
;

def subtract_(A, B):
    return [i for i in A if i not in B];
;

def reduce_octant(A):
    c = 0;
    while c < len(A):
        A = subtract_(A,
                subtract_([[A[c][0], -A[c][1]],
                 [-A[c][0], A[c][1]],
                 [-A[c][0], -A[c][1]],
                 [A[c][1], A[c][0]],
                 [A[c][1], -A[c][0]],
                 [-A[c][1], A[c][0]],
                 [-A[c][1], -A[c][0]],
                ], [A[c]])
            );
        c = c + 1;
    ;
    return A;
;

def all_distinct(L):
    flag = True;
    for i in L:
        if L.count(i) > 1:
            flag = False;
            break;
    return flag;
;

def find_triangle():
    G1 = direct_product_permgroups([SymmetricGroup(2),SymmetricGroup(2)]);
    G2 = SymmetricGroup(2);

    triangles = [];
    # the unit circle in (F_q)^2
    F2 = unit_circle();
    O = [0, 0];

    for P in F2:
        w = walltime();
        C = P;
        nextF2 = add_set(C, F2);
        for Q in nextF2:
            if quadrance(Q, O) == 1:
                triangles = triangles + [[O, P, Q]];
        print("vertex:", P, "time:", walltime(w));
    return triangles;
;

def find_equi_quad():
    quads = [];
    F2 = unit_circle();
    O = [0, 0];

    for P in F2:
        w = walltime();
        C = P;
        nextF2 = add_set(C, F2);
        for Q in nextF2:
            C = Q;
            nextF2 = add_set(C, F2);
            if O != Q:
                for R in nextF2:
                    if quadrance(R, O) == 1 and P != R and angle(O, P, Q) == angle(P, Q, R):
                        quads = quads + [[O, P, Q, R]];
        print("vertex:", P, "time:", walltime(w));
    return quads;
;

def find_equi_pent():
    pents = [];
    F2 = unit_circle();
    O = [0, 0];

    for P in F2:
        w = walltime();
        C = P;
        nextF2 = add_set(C, F2);
        for Q in nextF2:
            C = Q;
            nextF2 = add_set(C, F2);
            if O != Q:
                for R in nextF2:
                    if O != R and P != R:# and angle(O, P, Q) == angle(P, Q, R):
                        C = R;
                        nextF2 = add_set(C, F2);
                        for S in nextF2:
                            if all_distinct([O, P, Q, R, S]):
                                if quadrance(S, O) == 1:#  and angle(P, Q, R) == angle(Q, R, S) == angle(R, S, O) == angle(S, O, P):
                                    pents = pents + [[O, P, Q, R, S]];
        print("vertex:", P, "time:", walltime(w));
    return pents;
;

def find_triangle_3d():

    triangles = [];
    H3 = octant_sphere();
    F3 = unit_sphere(H3);
    O = [0, 0, 0];

    for P in F3:
        w = walltime();
        C = P;
        nextF3 = add_set(C, F3);
        for Q in nextF3:
            if dot(Q, Q) == 1:
                triangles = triangles + [[O, P, Q]];
        print("vertex:", P, "time:", walltime(w));
    return triangles;
;

def find_equi_quad_3d():
    quads = [];
    H3 = octant_sphere();
    F3 = unit_sphere(H3);
    O = [0, 0, 0];

    for P in F3:
        w = walltime();
        C = P;
        nextF3 = add_set(C, F3);
        for Q in nextF3:
            C = Q;
            nextF3 = add_set(C, F3);
            if O != Q:
                for R in nextF3:
                    quad = [O, P, Q, R];
                    if quadrance(R, O) == 1 and all_distinct(quad): # and angle(O, P, Q) == angle(P, Q, R):
                        quads = quads + [quad];
        print("vertex:", P, "time:", walltime(w));
    return quads;
;

def find_tetra():
    quads = [];
    H3 = octant_sphere();
    F3 = unit_sphere(H3);
    O = [0, 0, 0];

    for P in F3:
        w = walltime();
        C = P;
        nextF3 = add_set(C, F3);
        for Q in nextF3:
            C = Q;
            nextF3 = add_set(C, F3);
            if quadrance(O, Q) == 1:
                for R in nextF3:
                    quad = [O, P, Q, R];
                    if quadrance(R, O) == 1 and quadrance(R, P) == 1: # and angle(O, P, Q) == angle(P, Q, R):
                        quads = quads + [quad];
        print("vertex:", P, "time:", walltime(w));
    return quads;
;

def frob(P, n):
    Q = [];
    for i in range(len(P)):
        P_ = [];
        for j in range(len(P[0])):
            P_.append((P[i][j])^(p^n));
        Q.append(P_);
    return Q;

def sym_rfl(P, n): # reflection along axes/planes
    Q = [];
    for i in range(len(P)):
        P1 = [];
        for j in range(len(P[0])):
            P1.append((-1 if 2*(j+1) != n(2*(j+1)) else 1) * P[i][j])
        Q.append(P1);
    return Q;
;

def sym_perm(P, n): # permutation of coordinates
    Q = [];
    for i in range(len(P)):
        P1 = [];
        for j in range(len(P[0])):
            P1.append(P[i][n(j+1)-1])
        Q.append(P1);
    return Q;
;

def list_conj(T):
    conj = [];
    for h in range(0, F.order().factor()[0][1]):
        R = frob(T, h);
        for i in range(len(R)):
            P = [];
            for j in range(len(R)):
                P.append(sub(R[Mod(i+j, len(R))], R[i]));
            ;
            for g in G1:
                for h in G2:
                    c_ = sym_perm(sym_rfl(P, g), h);
                    conj = conj + [c_];
    ;
    return conj;
;

def sieve_rotoreflect(S):
    multi = [0] * 150;
    i = 0;
    S = sorted(S);
    while i < len(S):
        w = walltime();
        conj = list_conj(S[i]);
        count = conj.count(S[i]);
        multi[count] = multi[count] + 1;
        conj = subtract_(conj, [S[i]]); # self conjugate?
        S = subtract(S, conj);
        print("grain:", S[i], "time:", walltime(w), " remaining:", len(S) - i);
        i = i + 1;
    return [S, multi];
;

def find_in_sorted_list(L, elem):
    i = bisect.bisect_left(L, elem)
    if i != len(L) and L[i] == elem:
        return i;
    return -1;
;

def subtract(A, B):
    for i in B:
        e = find_in_sorted_list(A, i);
        if e >= 0:
            A.pop(e);
    return A;
;

def sparsify(array):
    res = [];
    for i in range(len(array)):
        if array[i] != 0:
            res.append([i, array[i]]);
    return res;
;

final = "";
final2 = "";
conj_count = "";

#G1 = direct_product_permgroups([SymmetricGroup(2),SymmetricGroup(2)]);
#G2 = SymmetricGroup(2);
G1 = direct_product_permgroups([SymmetricGroup(2),SymmetricGroup(2),SymmetricGroup(2)]);
G2 = SymmetricGroup(3);
#G1 = direct_product_permgroups([SymmetricGroup(2),SymmetricGroup(2),SymmetricGroup(2),SymmetricGroup(2)]);
#G2 = SymmetricGroup(4);

for p_ in range(3, 42):
    if is_prime(p_) == 1:
        v = walltime();
        p = p_;
        F = FiniteField(p);

        S = find_tetra();
        print("find time:", walltime(v));
        v = walltime();
        final = final + "(" + str(p_) + ", " + str(len(S)) + "), ";

        S = sieve_rotoreflect(S);

        print("sieve time:", walltime(v));
        final2 = final2 + "(" + str(p_) + ", " + str(len(S[0])) + "), ";

        conj_count = conj_count + "(" + str(p_) + ", " + str(sparsify(S[1])) + "), ";
        print("* * *");
;

print(final);
print(final2);
print(conj_count);
