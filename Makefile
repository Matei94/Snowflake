.PHONY: build
build: snowflake.c helpers.h helpers.c
	gcc -Wall snowflake.c helpers.c -lX11 -lm -o snowflake

.PHONY: run
run: snowflake
	./snowflake

.PHONY: clean
clean:
	rm snowflake snowflake.ppm
