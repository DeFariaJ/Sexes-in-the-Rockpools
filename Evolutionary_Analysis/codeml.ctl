seqfile = z.phy    * sequence data file name
     treefile = tree.tree  * tree structure file name
      outfile = mlc          * main result file name
        noisy = 0   * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
      seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree
       aaDist = 0   * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
        model = 0
      NSsites = 0 1 2 3 7 8
                    * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                    * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
        icode = 0   * 0:standard genetic code; 1:mammalian mt; 2-10:see below
        Mgene = 0   * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all
    fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2   * initial or fixed kappa
    fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate 
        omega = 1  * initial or fixed omega, for codons or codon-based AAs
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
    RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
    Small_Diff = .5e-6 * small value used in the difference approximation of derivatives
    cleandata = 0   * remove sites with ambiguity data (1:yes, 0:no)?
    fix_blength = 1   * 0: ignore, -1: random, 1: initial, 2: fixed
      method = 0   * 0: simultaneous; 1: one branch at a time