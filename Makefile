all:
	make -C ./cotocoa/src
	make -C ./cotocoa/test
	make -C ./main/coupler
	make -C ./main/requester
	make -C ./main/worker
	@echo "ビルド成功"
