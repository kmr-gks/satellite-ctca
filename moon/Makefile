.PHONY: all cotocoa coupler requester worker clean

all: cotocoa coupler worker requester
	@echo "ビルド成功"

cotocoa:
	make -C ./requester/lib/cotocoa/src

coupler:
	make -C ./coupler

requester:
	cotocoa_mode=COTOCOA_REQUESTER ./build.sh || { echo "Build failed in requester"; exit 1; }

worker:
	make -C ./worker

clean:
	make -C ./requester/lib/cotocoa/src clean
	make -C ./coupler clean
	make -C ./requester clean
	make -C ./worker clean
	@echo "クリーン成功"
