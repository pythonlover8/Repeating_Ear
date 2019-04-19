optimize=-o
FFLAGS=$(optimize) 
#CFLAGS=-g 
LDFLAGS=-shared
DEBUGFLAGS=-O0 -D _DEBUG
SACLIB=-L$(SACHOME)/lib -lsac -lsacio

Relo = bootstrap ellpse forward hist 
Relo1 = $(Relo) normalize stack takeoff xcorr_shift
Dep = matrix.c model_search.c nrutil.c weight.c sacio.c 

.PHONY: all
all: $(Relo1)

bootstrap: bootstrap.c 
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ $(Dep) \
	$(LDFLAGS) 
	@echo "$$< compiled" 

ellpse: ellpse.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o  $@ $(Dep) \
	$(LDFLAGS) 
	@echo "$$< compiled" 

forward: forward.c 
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o  $@ $(Dep) \
	$(LDFLAGS) 
	@echo "$$< compiled"

#grid_search: grid_search.c
#	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ \
#	$(LDFLAGS) 
#	$(ECHO) "$@ compiled" 

hist: hist.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ $(Dep) \
	$(LDFLAGS) 
	@echo "$$< compiled"  

normalize: normalize.c 
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ $(Dep)\
	$(LDFLAGS) -L$(SACHOME)/lib -lsac -lsacio 
	@echo "$$< compiled" 

stack: stack.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ $(Dep)\
	$(LDFLAGS) -L$(SACHOME)/lib -lsac -lsacio
	@echo "$$< compiled"

takeoff: takeoff.c 
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ $(Dep)\
	$(LDFLAGS) -L$(SACHOME)/lib -lsac -lsacio
	@echo "$$< compiled" 

xcorr_shift: xcorr_shift.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) $< -o $@ $(Dep)\
	$(LDFLAGS) -L$(SACHOME)/lib -lsac -lsacio
	@echo "$$< compiled"

.PHONY: clean
clean:
	rm -rf $(Relo1).o
