CC     = gcc
AR     = ar
CFLAGS = -g -Wall
# -std=gnu11 travis complains about this
# -std=c99 or -std=gnu99 might be important for compilation on a mac, otherwise duplicate definition for inline

SOURCES := $(wildcard **/*.c)
OBJECTS := $(patsubst %.c, %.o, $(SOURCES))

ifeq (1, $(CF_OPTIMIZE))
	CFLAGS += -O2
	CFLAGS := $(filter-out -g,$(CFLAGS))
endif

main: libgsl.a

%.o : %.c
	$(CC) -c $(CFLAGS) -I. $< -o $@

# main: $(OBJECTS)
# 	@echo $(OBJECTS)

libgsl.a: $(OBJECTS)
	@-rm -f $@
	$(AR) -csr $@ $^

test: libgsl.a
	gcc -I. test/test.c libgsl.a -std=c99 -lm -o test-main

clean:
	rm -f $(OBJECTS)

purge: clean
	rm -f libgsl.a

cleanhist:
	rm -rf .git .gitignore

