# Errata
As time strives by, some tragic errors were found by the author; while not having an horrific impact on the key conclusions, they're still extremly annoying,
and the author apologizes dearly for them. A short list of them is given below, which may (hopefully not) be expanded in the future:\

In Thesis:\
p. 26 the last integrand is sin(k*(r sin(phi)+r*cos(phi)))/(k*(sin(phi) cos(phi))), the division is missing; however, it was correctly implemented in the code.\
p. 28 first line: at the very end of the line a factor 2 is too many.\
p. 33 formula (2.71): A m-1 should be on top of the product sign instead of the m.\
p. 40 whole page: all F should be F^q.\
p. 77 second last line: also a factor 2 too many; instead of  [-1.3 pi, 1.3 pi], [-0.65 pi, 0.65 pi] is meant.\

In Code:\
File Chapter 7 - 4D Kuramoto with polynomial estimators: In line 26 a sign error occured. This was corrected. A new calculation yields for the 4D Kuramoto model for the
                                              indicator P3, that it's positive in 45% of all cases and in 7.4% of these cases it still converges to the sync state. Therefore
                                              P3 is not an "utter disappointment" anymore. In fact, it is possible to construct a quite good estimator for the sync basin
                                              by P1 || P2 || P3, which is positive in 73% of all caases and does still converge towwards  [0 0 0 0] in only 4-5% of those
