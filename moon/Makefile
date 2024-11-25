all:
	make -C ../cotocoa/src
	make -C ../cotocoa/test
	make -C ./coupler
	make -C ./requester
	make -C ./worker
	@echo "ビルド成功"

clean:
	make -C ../cotocoa/src clean
	make -C ../cotocoa/test clean
	make -C ./coupler clean
	make -C ./requester clean
	make -C ./worker clean
	@echo "クリーン成功"
