#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <arpa/inet.h>

void print_ipaddr()
{
    int fd;
    struct sockaddr_in dst_addr={0};
    struct sockaddr_in src_addr={0};
    socklen_t addrlen;
    char str[16];

    fd=socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP);
    dst_addr.sin_family=AF_INET;
    dst_addr.sin_port=htons(7);
    inet_aton("10.0.0.0",&dst_addr.sin_addr);
    connect(fd,(struct sockaddr *)&dst_addr,sizeof(dst_addr));
    addrlen=sizeof(src_addr);
    getsockname(fd,(struct sockaddr *)&src_addr,&addrlen);
    inet_ntop(AF_INET,&src_addr.sin_addr,str,sizeof(str));
    printf("%s\n",str);
    close(fd);
}

int sockfd;
int init_server()
{
    int client_sockfd;
    struct sockaddr_in addr;
 
    socklen_t len = sizeof( struct sockaddr_in );
    struct sockaddr_in from_addr;
 
    print_ipaddr();

    if( ( sockfd = socket( AF_INET, SOCK_STREAM, 0 ) ) < 0 ) {
        perror( "socket" );
    }
 
    addr.sin_family = AF_INET;
    addr.sin_port = htons( 10000);
    addr.sin_addr.s_addr = INADDR_ANY;

    int yes = 1;
    if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, (const char *)&yes, sizeof(yes)) < 0) {
        perror( "setsockopt" );
    }

    if( bind( sockfd, (struct sockaddr *)&addr, sizeof( addr ) ) < 0 ) {
        perror( "bind" );
    }
 
    if( listen( sockfd, SOMAXCONN ) < 0 ) {
        perror( "listen" );
    }
 
    if( ( client_sockfd = accept( sockfd, (struct sockaddr *)&from_addr, &len ) ) < 0 ) {
        perror( "accept" );
    }

    return client_sockfd;
}

int finish_server(int *client_sockfd)
{
    close( *client_sockfd );
    close( sockfd );
    return 0;
}

int send_real4(int *sock, float *dat, int *len)
{
  send(*sock, dat, *len*sizeof(float), 0);
  return *len;
}

int recv_real4(int *sock, float *dat, int *len)
{
  int rsize;

  rsize = recv(*sock, dat, *len*sizeof(float), 0);

  return rsize;
}
