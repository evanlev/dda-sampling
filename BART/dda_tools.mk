# Copyright 2014. The Regents of the University of California.
# All rights reserved. Use of this source code is governed by 
# a BSD-style license which can be found in the LICENSE file.


dda_toolssrcs := $(wildcard $(srcdir)/dda_tools/*.c)
dda_toolscudasrcs := $(wildcard $(srcdir)/dda_tools/*.cu)
dda_toolsobjs := $(dda_toolssrcs:.c=.o)

ifeq ($(CUDA),1)
dda_toolsobjs += $(dda_toolscudasrcs:.cu=.o)
endif


.INTERMEDIATE: $(dda_toolsobjs)

lib/libdda_tools.a: libdda_tools.a($(dda_toolsobjs))




