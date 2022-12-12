FROM ubuntu:22.04
RUN apt-get update --fix-missing \
&& apt-get install -y wget unzip libgomp1 \
&& cd /usr/local/ \
&& wget https://github.com/alexdobin/STAR/releases/download/2.7.10a_alpha_220818/STAR_2.7.10a_alpha_220818_Linux_x86_64_static.zip \
&& unzip STAR_2.7.10a_alpha_220818_Linux_x86_64_static.zip \
&& rm STAR_2.7.10a_alpha_220818_Linux_x86_64_static.zip \
&& cp STAR /usr/local/bin/ \
&& chmod +x /usr/local/bin/STAR \
&& rm STAR

ENTRYPOINT ["/usr/local/bin/STAR"]
