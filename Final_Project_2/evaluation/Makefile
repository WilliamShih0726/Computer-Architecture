make_distribute:
	g++ -g distribute.cpp -lpthread -o distribute
run_distribute:
	./distribute ../data/gencode.v43.transcripts.noN.fa ../data/queries.10K.fa 11 256 1
make_index:
	/opt/riscv-ca-fp2/bin/riscv64-unknown-elf-gcc -static -I /home/gem5/include -o index.o index.c /home/gem5/util/m5/build/riscv/out/libm5.a

run_index:
	/home/gem5/build/RISCV/gem5.opt /home/gem5/configs/example/se.py -c index.o -o "11" --cpu-type=MinorCPU --cpu-clock=800MHz --mem-type=Ramulator --ramulator-config=/home/gem5/configs/ramulator/LPDDR4-config.cfg --mem-size=256MB --caches --l1i_size=32kB --l1i_assoc=8 --l1d_size=32kB --l1d_assoc=8 --cacheline=32
	rm -rf index_result/m5out
	mv m5out index_result/

make_quantify:
	/opt/riscv-ca-fp2/bin/riscv64-unknown-elf-gcc -static -I /home/gem5/include -o quantify.o quantify.c /home/gem5/util/m5/build/riscv/out/libm5.a

run_quantify:
	/home/gem5/build/RISCV/gem5.opt /home/gem5/configs/example/se.py -c quantify.o -o "11" --cpu-type=MinorCPU --cpu-clock=800MHz --mem-type=Ramulator --ramulator-config=/home/gem5/configs/ramulator/LPDDR4-config.cfg --mem-size=256MB --caches --l1i_size=32kB --l1i_assoc=8 --l1d_size=32kB --l1d_assoc=8 --cacheline=32
	rm -rf quantify_result/m5out
	mv m5out quantify_result/

clean_distribute:
	rm distribute
	rm transcripts_data.h
	rm query_data.h
clean_index:
	rm index.o
	rm index_result/hash-map.txt
	rm -rf index_result/m5out/
clean_quantify:
	rm quantify.o
	rm quantify_result/final_result.txt
	rm -rf quantify_result/m5out/