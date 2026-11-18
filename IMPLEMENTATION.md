# Γ₁(N) support in `ArithGroupRQ` — implementation plan

## Goal

`ArithGroupRQ(F, level=N)` should model

    Γ₁(N) = { [[a, b], [c, d]] ∈ SL₂(O_F) : c ≡ 0, d ≡ 1 (mod N) }

for an integral ideal `N` of a real quadratic field `F` (`level=1` by default,
which must reproduce today's behaviour for SL₂(O_F) exactly). As for level 1, the
class must provide:

- a finite list of generators (`self.Ugens`, `self._gens`),
- a complete list of relators (`self._relation_words`), and
- a word-problem solver `get_word_rep(gamma)`.

All the darmonpoints machinery used downstream (`abelianization`,
`find_bounding_cycle`, `OneChains`, …) needs only those three.

## Approach: hand-rolled Reidemeister–Schreier + GAP for Tietze simplification

The current level-1 code already gives a presentation ⟨T₁, T_w, U, S | R⟩ of
G = SL₂(O_F) and a word-problem solver for it (continued fractions). Γ₁(N) has
finite index in G, and the right cosets Γ₁(N)\G are explicit:

    Γ₁(N)·g  ↦  bottom row of g mod N  ∈  { (c, d) ∈ (O/N)² : (c, d) + N = O }

because left multiplication by h ∈ Γ₁(N) does not change the bottom row of g mod N.
G acts on these cosets on the right by `v ↦ v·g` (row vector times matrix,
reduced mod N). So the permutation action of every generator on the cosets is
cheap to compute, and Reidemeister–Schreier can be run directly.

**Why not do it all in GAP?** In GAP we would have to define Γ₁(N) inside the fp
group either from subgroup generators (which we don't have, since finding them is
the whole point) or as the preimage of a point stabiliser under a homomorphism to
a permutation group, and then use `IsomorphismFpGroup`. That works, but then both
the coset table and the rewriting are opaque. With the hand-rolled version:

- the coset table comes straight from the arithmetic, not from coset enumeration;
- rewriting an ambient word into Schreier generators takes time linear in the word
  length, and it also tests membership (the walk has to return to the base coset);
- every step can be checked against matrices.

GAP (through Sage's `libgap`) is only used for what it is best at: Tietze
simplification of the resulting (large) presentation, keeping track of the
isomorphism so that words can be translated.

## Steps

### 1. Refactor the level-1 data into "ambient" data
- Keep the current generators, relators and word solver, under new names:
  `self._ambient_Ugens`, `self._ambient_relation_words`,
  `self._get_word_rep_ambient(gamma, check)` (the current `get_word_rep` body).
- `level` is trivial → `self.Ugens`, `self._gens`, `self._relation_words` are the
  ambient ones and `get_word_rep` delegates to the ambient solver. No behaviour
  change for existing scripts (`examples/main.sage`, `examples/test_carlos.sage`).

### 2. Normalise the level
- `self.level = F.ideal(level)`; require it to be integral and nonzero.
- Coset keys: `(N.reduce(c), N.reduce(d))` (canonical residues mod N).

### 3. Coset table and transversal (BFS)
- Start at v₀ = (0, 1) (the coset Γ₁(N) itself). Run a BFS over cosets using the
  generators and their inverses.
- Record:
  - `self._cosets`: list of keys (index 0 = v₀), plus a dict from key to index;
  - `self._coset_table[k][i]` = index of v_k·g_i, and `self._coset_table_inv[k][i]` = index of v_k·g_i⁻¹;
  - transversal words `t_k` (ambient Tietze words) and matrices `r_k`, with
    `(0, 1)·r_k ≡ v_k (mod N)`;
  - the tree edges: the pairs (k, i) with t_k·g_i = t_{k'} freely, meaning
    k' = table[k][i] was discovered from k through the letter +i (or the
    symmetric case for −i).
- Sanity check: the number of cosets equals
  `[G : Γ₁(N)] = Nm(N)² ∏_{P | N} (1 − Nm(P)⁻²)`.

### 4. Schreier generators
- For each pair (k, i) that is not a tree edge: s_{k,i} = r_k · g_i · r_{k'}⁻¹, where k' = table[k][i].
  Number them 1..m. There are m = [G:Γ₁(N)]·(#gens − 1) + 1 of them; this is checked.
- Check that every s_{k,i} lies in Γ₁(N).

### 5. Reidemeister rewriting `_rewrite_schreier(word, start=0)`
- Walk the word, keeping track of the current coset k:
  - letter +i: emit s_{k,i} (unless (k, i) is a tree edge), then k ← table[k][i];
  - letter −i: k' ← table_inv[k][i], emit s_{k',i}⁻¹ (unless (k', i) is a tree edge), then k ← k'.
- Return the emitted word and the final coset.
- Relators of Γ₁(N): `_rewrite_schreier(R, start=k)` for every ambient relator R
  and every coset k. Drop the empty ones.

### 6. Tietze simplification (GAP)
- Build `Q = FreeGroup(m) / relators` in libgap and call
  `iso = IsomorphismSimplifiedFpGroup(Q)`. Option `simplify=True` (default);
  `simplify=False` keeps the raw Schreier presentation, which is useful for
  debugging and as a fallback.
- From `iso`, extract:
  - for each new generator, a preimage word in the Schreier generators
    (`PreImagesRepresentative`), which gives its matrix in `self.Ugens`;
  - for each Schreier generator, its image word in the new generators
    (`ImageElm`), stored as `self._schreier_to_gens`;
  - the relators of `Range(iso)`, stored as `self._relation_words`.

### 7. Word problem for Γ₁(N)
`get_word_rep(gamma, check=True)`:
1. Membership test: entries integral, det = 1, c ∈ N, d − 1 ∈ N. Otherwise raise `ValueError`.
2. `w = self._get_word_rep_ambient(gamma)`.
3. `(s, k_end) = self._rewrite_schreier(w)`; assert `k_end == 0`.
4. Substitute each Schreier letter by its image under `self._schreier_to_gens`
   and freely reduce (`reduce_word_tietze`).
5. If `check`, multiply out and compare with gamma.

### 8. Consistency checks at construction time (as for level 1)
- Every relator evaluates to the identity matrix.
- Every generator lies in Γ₁(N).

## API summary

```python
G = ArithGroupRQ(F, level=N, simplify=True, use_shapiro=False)
G.level            # ideal N
G.index()          # [SL2(O_F) : Γ₁(N)]
G.transversal()    # transversal matrices r_k (not named coset_reps, a name that
                   # darmonpoints uses for S-arithmetic groups)
G.gens(), G.get_relation_words(), G.get_word_rep(gamma)
```

## Testing (`tests/test_arithgrouprq.py`, F = Q(√5))

Levels: (1), (2) (inert, index 15), (√5) (ramified, index 24),
(3/2·√5 − 1/2) (split prime above 11, index 120), and (2)² (index 240) if it runs fast enough.

1. Level 1 is unchanged: same generators, and the word for `gamma4` multiplies back.
2. The coset count matches the index formula.
3. All relators are trivial, and all generators lie in Γ₁(N).
4. Word problem round trip on many random elements of Γ₁(N). To get them, take a
   random γ ∈ SL₂(O_F) and replace it by γ · r_k⁻¹, where k is its coset.
5. An element outside Γ₁(N) raises `ValueError`.
6. Independent check, via GAP coset enumeration in the ambient fp group:
   `Index(G_fp, Subgroup(G_fp, new gens as ambient words))` equals the index.
   This shows the generators really generate Γ₁(N), relative to the ambient
   presentation.
7. The abelian invariants of the raw Schreier presentation and of the simplified
   one agree, and `G.abelianization()` runs.
8. Smoke test: `find_bounding_cycle` on an element of Γ₁(N), checked as in
   `examples/test_carlos.sage`.

## Finding: the original SL₂(O_F) presentation was incomplete

Test 6 failed with the original 11 ambient relators. Take the Γ₁(N) generators
and map them into the ambient fp group G_fp in two ways:

- through the Schreier words t_k·g_i·t_{k'}⁻¹: GAP's coset enumeration gives the
  correct index;
- through the continued-fraction word solver: GAP gives twice the index (30
  instead of 15 for N = (2), and similarly for (√5) and (3)).

Both routes give the same matrices, so G_fp has nontrivial elements that are the
identity in SL₂(O_F). In other words, the 11 relators present a proper cover of
SL₂(O_F).

The missing relator is the unit-ε analogue of (ST)³ = S², i.e. Cohn's universal
GE₂ relation for the unit ε:

    T_ε · S · T_{ε⁻¹} · S · T_ε · S = S² · U

It is now part of `_ambient_relation_words`. With it, test 6 passes for
N = (2), (√5), (3), 𝔭₁₁, and H₁ of the ambient fp group changes from ℤ/2 to 0.
This also changes level-1 results: `G.abelianization()` for SL₂(O_F) is now
trivial, so e.g. `min_exp` in `examples/main.sage` returns 1 instead of 2. Bounding
cycles computed before the change remain valid, since they only used true
relations.

Passing test 6 at these levels is strong evidence that the new ambient
presentation is complete, but it is not a proof.

## Results (F = Q(√5))

| N | index | Schreier gens | after Tietze | relators | H₁(Γ₁(N), ℤ) |
|---|---|---|---|---|---|
| (2) | 15 | 46 | 4 | 11 | (2, 2, 4, 4) |
| (√5) | 24 | 73 | 3 | 8 | (5, 5) |
| (3) | 80 | 241 | 4 | 21 | (3, 3, 3, 3) |
| 𝔭₁₁ = (3/2·√5 − 1/2) | 120 | 361 | 3 | 17 | (11, 11) |
| (4) | 240 | 721 | 4 | 31 | (4, 4, 8, 8) |

Construction takes a few seconds or less. Solving the word problem takes about
10–100 ms per element.

## Assumptions / caveats

- Correctness of the Γ₁(N) presentation depends on the ambient presentation of
  SL₂(O_F) being complete (see above). It is validated for Q(√5), as used in this
  repo; for other fields the ambient relators should be double-checked, for
  instance with test 6.
- The ambient word problem relies on quadratic continued fractions, so it needs
  O_F to be (norm-)Euclidean, as the current code already does.
- Raw presentations have about 3·[G:Γ₁(N)] generators. GAP's Tietze simplification
  handles indices in the low thousands without trouble; much larger levels may be slow.
