== colvarmodule.h ==

  /// \brief Number of metadynamics biases initialized (in normal
  /// conditions should be 1)
  static size_t n_meta_biases;
+  /// \brief Number of cvhd biases initialized (in normal
+  /// conditions should be 1)
+  static size_t n_cvhd_biases;

== colvarmodule.cpp ==

#include "colvarbias_meta.h"
+#include "colvarbias_cvhd.h"
#include "colvarbias_abf.h"

(...)

  /// initialize metadynamics instances
  parse_biases_type<colvarbias_meta>(conf, "metadynamics", n_meta_biases);
+
+  /// initialize cvhd instances
+  parse_biases_type<colvarbias_cvhd>(conf, "cvhd", n_cvhd_biases);

(...)

size_t                    colvarmodule::n_meta_biases = 0;
+size_t                    colvarmodule::n_cvhd_biases = 0;

== colvarbias.cpp ==

  if (to_lower_cppstr(key_str) == std::string("metadynamics")) {
    rank = cvm::n_meta_biases+1;
  }
+  if (to_lower_cppstr(key_str) == std::string("cvhd")) {
+    rank = cvm::n_cvhd_biases+1;
+  }


== Makefile.g++ ==

colvarcomp_rotations.o: colvarcomp_rotations.cpp colvarmodule.h \
 colvartypes.h colvarproxy.h colvarvalue.h colvarparse.h colvar.h \
 colvarcomp.h colvaratoms.h
+colvarcomp_bondbreak.o: colvarcomp_bondbreak.cpp colvarmodule.h \
+ colvartypes.h colvarproxy.h colvarparse.h colvarvalue.h colvaratoms.h \
+ colvar.h colvarcomp.h

(...)

colvarbias_meta.o: colvarbias_meta.cpp colvar.h colvarmodule.h \
 colvartypes.h colvarproxy.h colvarvalue.h colvarparse.h \
 colvarbias_meta.h colvarbias.h colvargrid.h
+colvarbias_cvhd.o: colvarbias_cvhd.cpp colvar.h colvarmodule.h \
+ colvartypes.h colvarproxy.h colvarvalue.h colvarparse.h \
+ colvarbias_cvhd.h colvarbias.h

== colvar.h  ==

  class selfcoordnum;
+  class bondbreak;
+  class angleswitch;

== colvar.cpp ==
  initialize_components("self-coordination "
                         "number",           "selfCoordNum",   selfcoordnum);
+
+  initialize_components("number of broken "
+                         "bonds",           "bondBreak",   bondbreak);
+
+  initialize_components("number of switched "
+                         "dihedrals",        "angleSwitch",   angleswitch);
