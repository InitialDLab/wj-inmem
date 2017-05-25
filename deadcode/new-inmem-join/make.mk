CXX=g++5
SRCDIR=../src
INCLUDEDIR=../include

SRC:=$(wildcard $(SRCDIR)/*.cpp)
INCLUDEFILES:=$(wildcard $(INCLUDEDIR)/*.h)

OBJS:=$(notdir $(SRC:.cpp=.o))
MAIN_OBJS:=$(filter main_%.o,$(OBJS))
OTHER_OBJS:=$(filter-out main_%.o,$(OBJS))
MAINS:=$(basename $(notdir $(MAIN_OBJS)))

INCLUDE=-I$(INCLUDEDIR)

.PHONY: all clean

all: $(MAINS)

$(MAINS): %: %.o $(OTHER_OBJS)
	$(CXX) $(CXXFLAG) $(INCLUDE) $^ -o $@ 

$(OBJS): %.o : $(SRCDIR)/%.cpp $(INCLUDEFILES)
	$(CXX) $(CXXFLAG) $(INCLUDE) -c $< -o $@

clean:
	rm -f $(MAINS) $(OBJS)



