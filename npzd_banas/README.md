# README for npzd_banas

#### This is edited npzd code to use with LiveOcean.

---

It is copied from the ROMS/Nonlinear/Biology folder of revision 1099, January 2022, after the introduction of varinfo.yaml.

I have also copied in ROMS/External/bio_Fennel.in and varinfo.dat.

Any code which I edit at all is first copied to ORIG_[].

The fundamental edits occur in:
- fennel.h where I make changes to the npzd functions to reproduce the Banas/Siedlecki/Davis model.  For each edited section I comment out the original code and then enclose the new code in ! PM Edit ... ! End PM Edit.
- bio_Fennel.in where I make changes to the parameter values.  We strictly retain the parameter names because they are baked into all the other code and so it is much easier to keep them.
- bio_Fennel_BLANK.in has a few more edits to the new bio_Fennel.in to allow it to be used in the dot_in code.  This mainly involved allowing for automated specification of boundary conditions, and turning on tracer sources and climatology.
- varinf0.yaml may need a few additional entries to work with the Fennel code.  I think this is a but in the ROMS source code.
