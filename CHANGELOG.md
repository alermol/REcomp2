# Change Log

## v0.0.1

---

## v2.0.0

* All codebase has been rewritten.
* Code become more parallel.

---

## v2.1.0

* Optimization of all steps.

---

## v2.1.1  

* Fix error with `style1.css`.
* Picture of graph layouts now is inside of html report so you can copy this file whenever you want and pictures will stay inside.

---

## v2.1.2

* `megablast` was replaced by `blastn` for removing junk 'other' contigs from fasta files.
* Progress bar was added on some steps instead less informative messages in standard output in terminal.
* Identity percent and query cover thresholds are not deteremenied by user now but calculated using clusterisation algorithms: K-means clustering (for runs with option `--include-other`) and agglomerative clustering (for runs without option `--include-other`).
* Dependence from `python_algorithms` package was eliminated.
* Dependeces from `tqdm` and `scikit-learn` has been added into `environment.yaml`
* Small changes in report.
* Minor fixes in code.
* Code refractoring.

---

## v2.1.2-1

* Test data become more lightweight

---

## v2.1.2-2

* Test data become more lightweight
* Small bugfixes

---

## v2.1.2-3

* Test data become more lightweight
* Small bugfixes

---

## v2.1.3

* Bug with cpu count was fixed
* Size of `REcomp2` was decreased
* `megablast` and `blastn` has been added as options for superclusters identification
* More accurate 'other' contigs selection algorithm based on artificial similarity index  
For each contig in superclusters will be calculated **Similarity Index** `(SI)`:

                                                    SI = abs(1 - (N / Q * 100))

    where N - number of hits for contig with all consensuses in supercluster, Q - maximal query cover per HSP across all HSP of contig. Contig with the smallest `SI` value will be selected.
