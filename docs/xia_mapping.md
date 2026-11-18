module	0	1	2	3	4	5	6	7	8
sid	2	3	4	5	6	7	8	9	10
type	500-14	500-14	100-14	100-14	100-14	100-14	100-14	100-14	100-14
number	1266	1226	1311	1313	1129	1286	1120	1115	1126
detector	TOF	T0D1-F	T0D1-B	T0D2-F0	T0D2-F1	Free	T0D2-B0	GAGG	T0D2-B1
ch0	t0d1-b15	t0d1-f0	t0d1-b0	t0d2-f0	t0d2-f16	t0d2-f22	t0d2-b0		t0d2-b16
ch1		t0d1-f1	t0d1-b1	t0d2-f1	t0d2-f17	t0d2-b1	error		t0d2-b17
ch2	Main_trigger	t0d1-f2	t0d1-b2	t0d2-f2	t0d2-f18	t0d2-b23	t0d2-b2	apd14	t0d2-b18
ch3		t0d1-f3	t0d1-b3	t0d2-f3	t0d2-f19	t0d2-b31	t0d2-b3	apd13	t0d2-b19
ch4		t0d1-f4	t0d1-b4	t0d2-f4	t0d2-f20		t0d2-b4	apd12	t0d2-b20
ch5		t0d1-f5	t0d1-b5	t0d2-f5	t0d2-f21		t0d2-b5	apd11	t0d2-b21
ch6		t0d1-f6	t0d1-b6	t0d2-f6	error		t0d2-b6		t0d2-b22
ch7		t0d1-f7	t0d1-b7	t0d2-f7	t0d2-f23		t0d2-b7	apd9	error
ch8		t0d1-f8	t0d1-b8	t0d2-f8	t0d2-f24		t0d2-b8	apd8	t0d2-b24
ch9		t0d1-f9	t0d1-b9	t0d2-f9	t0d2-f25		t0d2-b9	apd7	t0d2-b25
ch10		t0d1-f10	t0d1-b10	t0d2-f10	t0d2-f26		t0d2-b10	apd6	t0d2-b26
ch11		t0d1-f11	t0d1-b11	t0d2-f11	t0d2-f27		t0d2-b11	apd5	t0d2-b27
ch12		t0d1-f12	t0d1-b12	t0d2-f12	t0d2-f28		t0d2-b12	apd4	t0d2-b28
ch13		t0d1-f13	t0d1-b13	t0d2-f13	t0d2-f29		t0d2-b13	apd3	t0d2-b29
ch14		t0d1-f14	t0d1-b14	t0d2-f14	t0d2-f30		t0d2-b14	apd2	t0d2-b30
ch15	tof	tod1-f15	error	t0d2-f15	t0d2-f31		t0d2-b15	apd1	error




现在我需要将这个项目移植到其他的实验设置中：
- 首先暂时先不考虑所有关于vme的程序，只考虑XIA的decode，以及根据XIA的main_trigger进行事件重构
- 原程序使用了3个XIA机箱，目前我们只使用了一个机箱，并且mapping如上所显示，因此我们只需要读取一个机箱的数据，同时
重构配套的jinr_t0d1, jinr_t0d2, jinr_beam, jinr_gagg， jinr_trigger即可
- 在decode的过程中进行mapping，数据格式和之前保持一致，同时使用trigger重构事件的方法也和之前保持一致。
请提出你的理解，疑问和建议


新增：sid 2 ch3 为vme_trigger，需要像之前一样重构出jinr_t1_trigger.
理解：
1. 不要移除vme 的所有程序，不需要改动原样的保留即可，只改动我让你改动的

2.beam 暂时理解为 sid 2的 ch15

疑问：

1. beam 暂定为 sid 2的 ch15
2. 输出的 编号使用0-13，apd10不存在
3.是的
4. 沿用之前的处理方式
5. 结构不变，目前只读取main_trigger
6. 带 external_time

建议：

1. 不要使用一个程序，和之前的一样，每个探测器使用一个crush
2. 


疑问：
1. 和之前保持一致
2.beam 暂时理解为 sid 2的 ch15
3. 可以，暂时按照0处理
4. 各一版，程序中增加 -main -t1选项
5. 可以


注意事件对齐的sift 需要保留，之后我会引入vme的事件，目前保留相关接口



1. beam 目前只对应 sid2 ch15即可
2. 旧逻辑即可
3. 探测器暂时不要t1版本，我需要之后sift之后再给出t1版本


1. 直接改旧的 crush_t0d1/t0d2/beam/trigger 的映射
2. 对严格按照mapping即可
3. 结构和之前保持一致