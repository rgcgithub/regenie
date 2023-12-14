# Helper Classes From Remeta
The classes in this folder are used to store sparse compressed covariance/LD
matrices from SKAT for meta-analysis. Because they rely on HTSlib, which might
not be available on a user's system, they only get compiled if the WITH_HTSLIB
macro is set.