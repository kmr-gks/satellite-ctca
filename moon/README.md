CoToCoAを用いてEMSESのデータを共有するサンプルコード

requesterの内のコードを一部変更する必要がある。

`allcom.F90`
以下の変数を追加する

```fortran
    !共有するphiの配列
    real*8,allocatable    :: phi_data(:)
    !共有する配列のサイズ
    integer(kind=8)               :: phi_data_size=35
    !共有領域のエリアID
    integer               :: phi_areaid
```

`esses.F90`
メインループ部分(`esses_mainstep`を呼び出す前か後にこのコードを追加する)

```fortran
if(usecotocoa.eq.1) call send_data
```

`ictcar.F90`
ccinit内に追加
```fortran
!-------------------- 
    !共有する配列の領域を確保
    allocate(phi_data(phi_data_size))

!-------------------- 
    !エリアIDを取得
    call CTCAR_regarea_real8(phi_data,phi_data_size,phi_areaid)
```
次のサブルーチンを追加
```fortran
subroutine send_data
    use allcom
    implicit none
    
    !リクエストを送るときのデータ
    integer(kind=4) ::req_params(10)
    integer(kind=4) ::x,y,z,phi_shape(5),i
    integer ::from_rank=10

    if(myid.eq.from_rank) then
        phi_shape=shape(phi)
        x=phi_shape(2)
        y=phi_shape(3)
        z=phi_shape(4)
        phi_data = phi(1,1:phi_data_size,y/2,z/2,1)
        print*, "CTCArequester: phi_data=", phi_data
        !リクエスト時のデータを設定
        req_params(1)=from_rank
        req_params(2)=phi_data_size
        !リクエスト時にデータを送ることができる
        call CTCAR_sendreq(req_params,size(req_params))
    end if
end subroutine send_data
```



moonフォルダで`make`する

default-condition,dshield1,exp_surfaceのうちの一つ
シミュレーションフォルダに移動し、sbatchする


