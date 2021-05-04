
if (options_.nucleotide_constraints && i2r[R_nuc] != NucConst_[j]) {
    continue;
}

/* skip if we're close to the ends and j and i are too close */
if (options_.DEPflg && j - i == 1 && i <= nuclen_ - 1 && 
    Dep1_[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
    continue;
} 
if (options_.DEPflg && j - i == 2 && i <= nuclen_ - 2 && 
    Dep2_[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
    continue;
}

/* ======================================= */

if (options_.nucleotide_constraints == 1 && i2r[R2_nuc] != NucConst_[j - 1]) {
    continue;
}

if (options_.DEPflg && Dep1_[ii2r[L_nuc * 10 + L2_nuc]][i] == 0) {
    continue;
}
if (options_.DEPflg && Dep1_[ii2r[R2_nuc * 10 + R_nuc]][j - 1] == 0) {
    continue;
}

/* ======================================= */


if (options_.nucleotide_constraints && i2r[Lp_nuc] != NucConst_[p]) {
    continue;
}

if (options_.DEPflg && p == i + 1 && Dep1_[ii2r[L_nuc * 10 + Lp_nuc]][i] == 0) {
    continue;
}
if (options_.DEPflg && p == i + 2 && Dep2_[ii2r[L_nuc * 10 + Lp_nuc]][i] == 0) {
    continue;
}


/* ======================================= */


   if (options_.nucleotide_constraints && i2r[R_nuc] != NucConst_[j]) {
       continue;
   }
  
   /* skip if we're close to the ends and j and i are too close */
   if (options_.DEPflg && j - i == 1 && i <= nuclen_ - 1 && 
       Dep1_[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
       continue;
   } 
   if (options_.DEPflg && j - i == 2 && i <= nuclen_ - 2 && 
       Dep2_[ii2r[L_nuc * 10 + R_nuc]][i] == 0) {
       continue;
   }

/* ======================================= */


   if (options_.nucleotide_constraints && i2r[L2_nuc] != NucConst_[i + 1]) {
           continue;
       }

       if (options_.DEPflg && Dep1_[ii2r[L_nuc * 10 + L2_nuc]][i] == 0) {
           continue;
       }

/* ======================================= */

  if (options_.nucleotide_constraints && i2r[R2_nuc] != NucConst_[j - 1]) {
      continue;
  }

  if (options_.DEPflg && Dep1_[ii2r[R2_nuc * 10 + R_nuc]][j - 1] == 0) {
      continue;
  }
/* ======================================= */

  if (options_.nucleotide_constraints && i2r[Lp2_nuc] != NucConst_[p - 1]) {
      continue;
  }

  if (options_.DEPflg && Dep1_[ii2r[Lp2_nuc * 10 + Lp_nuc]][p - 1] == 0) {
      continue;
  }

/* ======================================= */

 if (options_.nucleotide_constraints && i2r[Rj1_nuc] != NucConst_[j - 1]) {
     continue;
 }

 if (options_.DEPflg && Dep1_[ii2r[Rj1_nuc * 10 + R_nuc]][j - 1] == 0) {
     continue;
 }
/* ======================================= */

  if (options_.nucleotide_constraints && i2r[Li1_nuc] != NucConst_[i + 1]) {
      continue;
  }
  if (options_.DEPflg && Dep1_[ii2r[L_nuc * 10 + Li1_nuc]][i] == 0) {
      continue;
  }
/* ======================================= */

 if (options_.nucleotide_constraints && i2r[Rj1_nuc] != NucConst_[j - 1]) {
     continue;
 }
 if (options_.DEPflg && Dep1_[ii2r[Rj1_nuc * 10 + R_nuc]][j - 1] == 0) {
     continue;
 }
/* ======================================= */

 if (options_.nucleotide_constraints && i2r[Lk_nuc] != NucConst_[k]) {
        continue;
    }
    if (options_.DEPflg && Dep1_[ii2r[Rk1_nuc * 10 + Lk_nuc]][k - 1] == 0) {
        continue;
