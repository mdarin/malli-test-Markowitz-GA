PROGNAME = alloc

$(PROGNAME): clean deps
	@echo "Building..."
	@./build.sh
	@echo "Done!"

deps:
	@echo "Gethering dependencies..."
	@./deps.sh
	@echo "Done!" 

clean:
	@echo "Cleaning..."
	@./clean.sh
	@echo "Done!"

clean_all: clean
	@echo "Sweepping up..."
	@./sweepup.sh
	@echo "Done!"

