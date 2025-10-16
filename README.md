# Background

This set of scripts began as an attempt to reverse-engineer a program called Color Vision Simulator developed for altering photographs to mimic various forms of primate color vision. The method used by CVS is described in the paper where it was introduced (Melin et al. 2013), but the software itself has not been published. I later added other functions for analyzing animal color vision more generally, and the image processing code is now contained in its own file (image_processing.py) because it relies on OpenCV while the main file uses matplotlib -- see [this bug report](https://github.com/opencv/opencv-python/issues/386).

I have no formal background in biology or neuroscience, so this project is meant for personal use. Professional-quality software with some similar functions includes [micaToolbox](https://www.empiricalimaging.com/) and [pavo](https://pavo.colrverse.com/) (both free to download).

# References

* A.D. Melin, D.W. Kline, C.M. Hickey, L.M. Fedigan, Food search through the eyes of a monkey: A functional substitution approach for assessing the ecology of primate color vision, Vision Research, Volume 86, 2013, Pages 87-96, ISSN 0042-6989, https://doi.org/10.1016/j.visres.2013.04.013. (https://www.sciencedirect.com/science/article/pii/S0042698913001119)
* Hogan JD, Fedigan LM, Hiramatsu C, Kawamura S, Melin AD. Trichromatic perception of flower colour improves resource detection among New World monkeys. Sci Rep. 2018 Jul 18;8(1):10883. doi: 10.1038/s41598-018-28997-4. PMID: 30022096; PMCID: PMC6052032. https://www.nature.com/articles/s41598-018-28997-4
* Melin AD, Chiou KL, Walco ER, Bergstrom ML, Kawamura S, Fedigan LM. Trichromacy increases fruit intake rates of wild capuchins (Cebus capucinus imitator). Proc Natl Acad Sci U S A. 2017 Sep 26;114(39):10402-10407. doi: 10.1073/pnas.1705957114. Epub 2017 Sep 11. PMID: 28894009; PMCID: PMC5625910. https://www.pnas.org/doi/10.1073/pnas.1705957114
* Fedigan LM, Melin AD, Addicott JF, Kawamura S. The heterozygote superiority hypothesis for polymorphic color vision is not supported by long-term fitness data from wild neotropical monkeys. PLoS One. 2014 Jan 3;9(1):e84872. doi: 10.1371/journal.pone.0084872. PMID: 24404195; PMCID: PMC3880319.
* Melin AD, Kline DW, Hiramatsu C, Caro T (2016) Zebra Stripes through the Eyes of Their Predators, Zebras, and Humans. PLoS ONE 11(1): e0145679. https://doi.org/10.1371/journal.pone.0145679
* Corso J, Bowler M, Heymann EW, Roos C, Mundy NI. Highly polymorphic colour vision in a New World monkey with red facial skin, the bald uakari (Cacajao calvus). Proc Biol Sci. 2016 Apr 13;283(1828):20160067. doi: 10.1098/rspb.2016.0067. PMID: 27053753; PMCID: PMC4843651. https://doi.org/10.1098/rspb.2016.0067
* Veilleux CC, Scarry CJ, Di Fiore A, Kirk EC, Bolnick DA, Lewis RJ. Group benefit associated with polymorphic trichromacy in a Malagasy primate (Propithecus verreauxi). Sci Rep. 2016 Dec 2;6:38418. doi: 10.1038/srep38418. PMID: 27910919; PMCID: PMC5133583. https://www.nature.com/articles/srep38418
* Hiramatsu C, Melin AD, Allen WL, Dubuc C, Higham JP. Experimental evidence that primate trichromacy is well suited for detecting primate social colour signals. Proc Biol Sci. 2017 Jun 14;284(1856):20162458. doi: 10.1098/rspb.2016.2458. PMID: 28615496; PMCID: PMC5474062. https://doi.org/10.1098/rspb.2016.2458
* Stavenga, D.G. On visual pigment templates and the spectral shape of invertebrate rhodopsins and metarhodopsins. J Comp Physiol A 196, 869–878 (2010). https://doi.org/10.1007/s00359-010-0568-7
* Hempel de Ibarra N, Vorobyev M, Menzel R. Mechanisms, functions and ecology of colour vision in the honeybee. J Comp Physiol A Neuroethol Sens Neural Behav Physiol. 2014 Jun;200(6):411-33. doi: 10.1007/s00359-014-0915-1. Epub 2014 May 15. PMID: 24828676; PMCID: PMC4035557.
* CMF From Cone Fundamentals. Horizon Lab @ UCRS. https://horizon-lab.org/colorvis/cone2cmf.html
* Petroc Sumner, Catherine A. Arrese, Julian C. Partridge; The ecology of visual pigment tuning in an Australian marsupial: the honey possum Tarsipes rostratus. J Exp Biol 15 May 2005; 208 (10): 1803–1815. doi: https://doi.org/10.1242/jeb.01610
* Arrese, A. C., Beazley, L. D. and Neumeyer, C. Behavioural evidence for marsupial trichromacy. doi:10.1016/j.cub.2006.02.036
* Ortín-Martínez A, Nadal-Nicolás FM, Jiménez-López M, Alburquerque-Béjar JJ, Nieto-López L, García-Ayuso D, Villegas-Pérez MP, Vidal-Sanz M, Agudo-Barriuso M. Number and distribution of mouse retinal cone photoreceptors: differences between an albino (Swiss) and a pigmented (C57/BL6) strain. PLoS One. 2014 Jul 16;9(7):e102392. doi: 10.1371/journal.pone.0102392. PMID: 25029531; PMCID: PMC4100816.
* Jacobs, G. H., Fenwick, J. A. and Williams, G. A. Cone-based vision of rats for ultraviolet and visible lights. Journal of Experimental Biology, 204(14), 15 July 2001. https://doi.org/10.1242/jeb.204.14.2439
* CIE 2-deg CMFs. http://cvrl.ucl.ac.uk/database/text/cienewxyz/cie2012xyz2.htm
* Pridmore, R.W. A new transformation of cone responses to opponent color responses. Atten Percept Psychophys 83, 1797–1803 (2021). https://doi.org/10.3758/s13414-020-02216-7
* https://chittkalab.sbcs.qmul.ac.uk/1992/Chittka%201992%20J%20Comp%20Physiol_R.pdf
* https://doi.org/10.1016/j.jqsrt.2020.107162
* https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12439

# Further reading
* https://xkcd.com/1926/
* http://xkcd.com/451/
