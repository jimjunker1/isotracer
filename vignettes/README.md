# Some run times

- Precompiling the vignettes on a 4-core laptop at 1700 MHz takes about 30 min.
- Precompiling the vignettes on a multicore desktop computer at 3600 MHz takes about 9 min.

# 2021-08-25

- Most of the vignettes in this folder are quite long to run, and several run processes in parallel for efficiency. This can lead to issues when running `R CMD check`.
- One solution, as suggested by a [post](https://blog.r-hub.io/2020/06/03/vignettes/) by Maëlle Salmon, is to pre-compute vignettes. Maëlle's post links to Jeroen Ooms's [post](https://ropensci.org/blog/2019/12/08/precompute-vignettes/) with detailed explanations about how to do it.
- As noted by Maëlle, this raises the question of testing and reproducibility. She cites [Henrik Bengtsson](https://www.mail-archive.com/r-package-devel@r-project.org/msg00812.html) who argues that testing can still be good if tests in `tests/` are good. For reproducibility, Henrik suggests to include the source vignettes in the package in `inst/full-vignettes/`.

