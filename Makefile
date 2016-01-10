build: snowflake.c
	gcc snowflake.c -lX11 -lm -o snowflake

run: snowflake
	./snowflake

clean:
	rm snowflake snowflake.ppm
