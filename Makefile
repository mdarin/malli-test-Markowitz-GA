PROGNAME = alloc

$(PROGNAME): clean_all deps
	@echo "Building..."
	@./bin/build.sh
	@echo "Done!"

deps:
	@echo "Gethering dependencies..."
	@./bin/deps.sh
	@echo "Done!" 

clean:
	@echo "Cleaning..."
	@./bin/clean.sh
	@echo "Done!"

clean_all: clean
	@echo "Sweepping up..."
	@./bin/sweepup.sh
	@echo "Done!"

