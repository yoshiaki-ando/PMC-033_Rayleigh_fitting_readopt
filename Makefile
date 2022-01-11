OBJS = Across_pts.o Geocoordinate.o Mie_scattering.o Msis.o across_point_atmosphere.o \
	associated_legendre.o beta.o calculate_intensity.o fitting_rayleigh.o \
	geometry.o get_observation_data.o get_rayleigh.o hcoord2geo.o interpolated_m.o \
	limb_point.o main.o optical_depth.o output_pmc.o pmc.o rayleigh_part.o \
	read_pmc_param.o search_pmc_region.o set_fitting_latlon.o solar_direction.o

HEADERS = Across_pts.h Date.h Geocoordinate.h Geoparameter.h Mathparameter.h Mie_scattering.h Msis.h Region.h interpolated_m.h pmc.h pmc_simulation.h

PREFIX = $(HOME)
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include

LIB = -lcomplex_bessel -lAndoLab -lnrlmsise_18
LINKER_OPTS = -Wl,-R$(LIB_DIR) -L$(LIB_DIR)

OPTS = -Wall -O3 -I$(INC_DIR)

.PHONY: all clean

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS) $(LINKER_OPTS) $(LIB)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

clean:
	rm -rf *.o main
