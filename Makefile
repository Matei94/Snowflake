.PHONY: build
build: snowflake.c
	gcc snowflake.c -lX11 -lm -o snowflake

.PHONY: run
run: snowflake
	./snowflake

.PHONY: clean
clean:
	rm snowflake snowflake.ppm
