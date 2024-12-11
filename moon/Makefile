all:
	cp -r ./requester-custom/* ./requester/src/cotocoa
	make -C ./requester/lib/cotocoa/src
	make -C ./coupler
	make -C ./requester requester
	make -C ./worker
	@echo "ビルド成功"

clean:
	make -C ./requester/lib/cotocoa/src clean
	make -C ./coupler clean
	make -C ./requester clean
	make -C ./worker clean
	@echo "クリーン成功"
