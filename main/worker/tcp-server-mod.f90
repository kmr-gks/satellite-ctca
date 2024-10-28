module tcp_server_mod
  interface
     function init_server() bind (c)
       integer :: init_server
     end function init_server

     function finish_server(sock) bind (c)
       integer :: sock
       integer :: finish_server
     end function finish_server

     function send_real4(sock, dat, len) bind(c)
       integer :: sock, len
       real(kind=4), dimension(1:len) :: dat
       integer :: send_real4
     end function send_real4

     function recv_real4(sock, dat, len) bind(c)
       integer :: sock, len
       real(kind=4), dimension(1:len) :: dat
       integer :: recv_real4
     end function recv_real4
  end interface

  public init_server, finish_server, send_real4, recv_real4
end module tcp_server_mod
