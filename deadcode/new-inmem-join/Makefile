.PHONY: all debug release clean cleanall

CXXFLAG=-Wall -std=c++1y

all: debug release

debug:
	@mkdir -p debug
	@cp make.mk debug/Makefile
	make -C debug CXXFLAG="-g $(CXXFLAG)"

release:
	@mkdir -p release
	@cp make.mk release/Makefile
	make -C release CXXFLAG="-O2 $(CXXFLAG)"

clean:
	make -C debug clean
	make -C release clean

cleanall:
	rm -rf debug release

