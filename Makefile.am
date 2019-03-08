# Name of executables: The prefix "noinst" means the executables don't
# have to be installed anywhere.
noinst_PROGRAMS =  circular_disc_rolling_clamp_ks circular_disc_full_clamp_ks
include $(top_srcdir)/config/makefile_templates/demo_drivers
#---------------------------------------------------------------------------

# Local sources that circular_cap_fvk code depends on:
circular_disc_rolling_clamp_ks_SOURCES = \
circular_disc_uniform_pressure_koiter_steigmann_rolling_clamp.cc

circular_disc_full_clamp_ks_SOURCES = \
circular_disc_uniform_pressure_koiter_steigmann_fully_clamped_out_pinned_in.cc

#---------------------------------------------------------------------------

circular_disc_rolling_clamp_ks_LDADD = -L@libdir@ -lC1_large_displacement_plate_models \
 -lC1_basis -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

circular_disc_full_clamp_ks_LDADD = -L@libdir@ -lC1_large_displacement_plate_models \
 -lC1_basis -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Local sources that Jack's own code depends on: This code also uses
# objects from Jack's own library. The source code for this library
# lives in user_src/jack_cool and gets compiled in to the jack_cool 
# library. If we use it, we have to link against it...
#jacks_own_code_SOURCES = jacks_own_code.cc

# Required libraries: Jack's own library ("jack_cool") and oomph-lib's 
# and "generic" and "poisson" libraries. Following the logic
# explained above, Jack's own library should go first since it may
# (or may not) depend on objects in the other libraries, whereas they will
# certainly not depend on objects in his. 
# jacks_own_code_LDADD = -L@libdir@ -ljack_cool -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#---------------------------------------------------------------------------


# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
# Automake will replace the variable @includedir@ with the actual
# include directory.

AM_CPPFLAGS += -I@includedir@ -I$(pathtoeigen) -I$(pathtospectra) -Wall #-g -DRANGE_CHECK -DPARANOID
#EXTRA_DIST = Circle1.1.ele Circle1.1.node Circle1.1.poly