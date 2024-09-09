import dcp
import msm
import zvp_glv_sac as zvpmsm
import utils
from sage.all import ZZ, EllipticCurve, GF, ceil, randint
from pari_tools.dcp_pari import glvdcpsolver_pari, dcpsolver_pari
from experiments import ZVPparams
from dcp import prepare_Vpolynomial, prepare_simpleVpolynomial


TEST_REGISTERS = utils.Registers()
X1, Y1, X2, Y2 = TEST_REGISTERS.gens
TEST_REGISTERS.add_tuple(X1 + X2 + 2, X1 + X2 + 2)


def test_glvdcpsolver_pari():
    secp256k1 = utils.GLVCurve()
    secp256k1.set_secp256k1()

    P = secp256k1.curve(
        102673692605890747955827626789425121715741455357436816214783425309028914628708,
        87043485697416879999829492078402404425907071118621005692711309620797730731601,
    )
    Q = secp256k1.curve(
        69160547611380020815779950572212268692227625588734315405709809992685171713095,
        49653363091173826186043652289215864926772053949623812607939578824633542393314,
    )
    k1, k2 = 3, 2
    registers = utils.Registers()
    X1, _, X2, _ = registers.gens
    decoy = X1 + X2
    polynomial = (
        X1**2
        + X2
        - 73228285294405723246573492566392993825226941833864860421887399373708380041916
    )
    registers.add_tuple(decoy, decoy)
    registers.add_tuple(polynomial, polynomial)

    Vs = [
        prepare_Vpolynomial(polynomial, registers, secp256k1)
        for polynomial in registers.polynomials
    ]
    assert Vs[1] == Vs[1].parent()(
        "3623740642296190556596022661315200510509771318140512627724956322491646257792400725600189268559480203350328236110123415136415655083604431178924234041790464*x^4*d1^2*n2^2 - 120395027177972605108971164049770150217768064901904679635358144052332456178816*x^4*n1*d1*n2*d2 + x^4*n1^2*d2^2 + 7247481284592381113192045322630401021019542636281025255449912644983292515584801451200378537118960406700656472220246830272831310167208862357848468083580928*x^2*n1*d1*n2^2 - 530720627173996962291869849190856726836471755115534911383876145372620410492084384246309376279373854583628241936423146376344945713842175756230048908653413237334904006900592346745909869758323382593774531669810287514581589968018178048*x^2*d1^2*n2^2 + 120395027177972605108971164049770150217768064901904679635358144052332456178816*x^2*n1^2*n2*d2 + 17632642796432617401468186556033965379842497992066244078915067620766194908190828026121813978651104228272071615700103407886635823043261051375769822542502912*x^2*n1*d1*n2*d2 - 146456570588811446493146985132785987650453883667729720843774798747416760083832*x^2*n1^2*d2^2 + 28*x^2*d1^2*d2^2 + 3623740642296190556596022661315200510509771318140512627724956322491646257792400725600189268559480203350328236110123415136415655083604431178924234041790464*n1^2*n2^2 - 530720627173996962291869849190856726836471755115534911383876145372620410492084384246309376279373854583628241936423146376344945713842175756230048908653413237334904006900592346745909869758323382593774531669810287514581589968018178048*n1*d1*n2^2 + 19431880749161692114439479801057861629732172516303529473218388436658014353212852273793871786640481279219300583231179012044401274084763740897259334334829173831724504553088943254015079590356681850630611504348168720624340706474169538803981060979833875411632726396411143778212986111277988811572102984546095529984*d1^2*n2^2 - 8816321398216308700734093278016982689921248996033122039457533810383097454095414013060906989325552114136035807850051703943317911521630525687884911271251456*n1^2*n2*d2 - 645604098595757822801450323956043684472783879251833997496610119181956055431013857684866811645551742935727475200395656289498595919889430303459892793967136039785348674874868904504529623046173993675056702518741179159450135247056029696*n1*d1*n2*d2 - 1685530380491616471525596296696782103048752908626665514895014016732654386503424*d1^2*n2*d2 + 5362381767158877501737147854532438911760514307467365184716554607310616836281697443898948406277520253223656572319414677181072963952857746521120913916951056*n1^2*d2^2 - 28*n1*d1*d2^2 - 2050391988243360250904057791859003827106354371348216091812847182463834641173648*d1^2*d2^2"
    )
    roots = glvdcpsolver_pari(
        secp256k1.p,
        secp256k1.curve.a4(),
        secp256k1.curve.a6(),
        k1,
        secp256k1.lam,
        k2,
        Vs,
        registers,
    )
    roots = [i[0] for i in roots]
    assert len(roots) != 0
    assert set(roots).issubset(
        set(
            [
                102673692605890747955827626789425121715741455357436816214783425309028914628708,
                86124371540462867129969031022427479710818463712531960032840963285688824050049,
            ]
        )
    ), set(roots)


def test_dcpsolver_pari():
    secp256k1 = utils.GLVCurve()
    secp256k1.set_secp256k1()
    beta = 60197513588986302554485582024885075108884032450952339817679072026166228089408
    lam = 78074008874160198520644763525212887401909906723592317393988542598630163514318
    P = secp256k1.curve(
        102673692605890747955827626789425121715741455357436816214783425309028914628708,
        87043485697416879999829492078402404425907071118621005692711309620797730731601,
    )
    k = 5
    registers = utils.Registers()
    X1, _, X2, _ = registers.gens
    decoy = X1 + X2
    polynomial = (
        X1**2
        + X2
        - 23507235017084783189512012573502341255999837337079867586857010079151143785546
    )
    registers.add_tuple(decoy, decoy)
    registers.add_tuple(polynomial, polynomial)
    Vs = [
        prepare_simpleVpolynomial(polynomial, secp256k1, registers, beta=1)
        for polynomial in registers.polynomials
    ]
    roots = dcpsolver_pari(secp256k1.field.order(), 0, 7, k, 1, Vs, registers)
    roots = [i[0] for i in roots]
    assert len(roots) != 0
    assert set(roots).issubset(
        set(
            [
                102673692605890747955827626789425121715741455357436816214783425309028914628708,
                55513075508145860462041760647794016639151178866915137189061803695425235878413,
            ]
        )
    ), set(roots)

    registers = utils.Registers()
    X1, _, X2, _ = registers.gens
    polynomial = (
        X1**2
        + X2
        - 40858653405389928731140344952206327443574865814776107614481074951822958796001
    )
    registers.add_tuple(decoy, decoy)
    registers.add_tuple(polynomial, polynomial)
    Vs = [
        prepare_simpleVpolynomial(polynomial, secp256k1, registers, beta=beta)
        for polynomial in registers.polynomials
    ]
    roots = dcpsolver_pari(secp256k1.field.order(), 0, 7, k, lam, Vs, registers)
    roots = [i[0] for i in roots]
    assert len(roots) != 0
    assert set(roots).issubset(
        set(
            [
                102673692605890747955827626789425121715741455357436816214783425309028914628708,
                99075450709825848198932375595499468931178626324563279679255491735703931462704,
            ]
        )
    ), set(roots)


def test_glvdcpsolver_pari2():
    secrets = utils.GLVSecrets()
    secrets.from_dict(
        {
            "k0": 9354505673776504235,
            "k1": 16488824988860969087,
            "k": 29909316888294194379575008222298337520,
            "curve": [0, 5],
            "field": 219367660988620325471586898010713677223,
            "beta": 46913353786251974458108963163358176822,
            "lambda": 92096968189742899780809622067096840422,
            "order": 219367660988620325450527123821457422649,
        }
    )
    glv = secrets.glv
    registers = utils.Registers()
    X1, Y1, X2, Y2, A, B = registers.all_gens
    f = X1**3 * X2**3 + X1**3 * B + X2**3 * B + B**2 - 1
    registers.add_tuple(*(f, Y1 * Y2 + 1))
    k = 2
    Vs = [
        prepare_simpleVpolynomial(polynomial, glv, registers, beta=1)
        for polynomial in registers.polynomials
    ]
    roots = dcpsolver_pari(
        glv.p, glv.curve.a4(), glv.curve.a6(), k, glv.lam, Vs, registers
    )
    # 65705721145317289516028210688079696622, 163286117200082233085479988588118356993
    P = glv.curve(
        65705721145317289516028210688079696622, 163286117200082233085479988588118356993
    )
    assert glv.curve(*roots[0]) == P
    Q = 2 * P
    assert P[1] * Q[1] + 1 == 0
    assert f(*P[:2], *Q[:2], glv.curve.a4(), glv.curve.a6()) == 0

    k1, k2 = -2, -2
    Vs = [
        prepare_Vpolynomial(polynomial, registers, glv)
        for polynomial in registers.polynomials
    ]
    roots = glvdcpsolver_pari(
        glv.p, glv.curve.a4(), glv.curve.a6(), k1, glv.lam, k2, Vs, registers
    )

    P1, Q1 = -P, (2 + 2 * glv.lam) * P
    assert (-2 - 2 * glv.lam) * P1 == Q1
    assert P1[1] * Q1[1] + 1 == 0
    assert f(*P1[:2], *Q1[:2], glv.curve.a4(), glv.curve.a6()) == 0
    R1, R2 = -2 * P1, -2 * P1
    assert Vs[0](P1[0], R1[0], 1, R2[0], 1) == 0


def test_solve_glv_dcp_pari():
    k1, k2 = ZZ(6), ZZ(6)
    glv = utils.GLVCurve()
    p = ZZ(1039)
    F = GF(p)
    glv.curve = EllipticCurve(F, [0, 6])
    P = glv.curve(47, 353)
    Q = 143 * P
    glv.order = ZZ(1033)
    glv.lam = ZZ(195)
    glv.beta = ZZ(140)
    assert P[0] + Q[0] + 2 == 0
    result = dcp.solve_glv_dcp_pari([k1, k2], glv, TEST_REGISTERS)
    assert int(result[0]) in [594, 47]


def test_solve_multi_dcp_pari():
    k, l = ZZ(-2), ZZ(3)
    glv = utils.GLVCurve()
    p = ZZ(664472600456690912117742264307)
    F = GF(p)
    glv.curve = EllipticCurve(F, [0, 7])
    P = glv.curve(422308990711158209733079302631, 148148080066598589241454378244)
    Q = glv.curve(93588305526752155239403665403, 111249103853384717460587115118)
    assert (l * P)[0] + Q[0] == 0
    registers = utils.Registers()
    X1, Y1, X2, Y2 = registers.gens
    registers.add(X1 + X2)
    result = dcp.solve_multi_dcp(k, l, glv, registers)
    assert int(result[0]) in [
        422308990711158209733079302631,
        168593139093414804842851950104,
        73570470652117897541811011572,
    ]


def test_solve_multi_dcp_pari2():
    k, l = ZZ(2), ZZ(-3)
    glv = utils.GLVCurve()
    p = ZZ(1203762391362186422796891378463)
    F = GF(p)
    glv.curve = EllipticCurve(F, [0, 7])
    P = glv.curve(744737971283560328960748838057, 654646481094124905305739376702)
    Q = glv.curve(521465934801467274552161980277, 1071698956917683922578329251129)
    assert (l * P)[0] + Q[0] == 0
    registers = utils.Registers()
    X1, Y1, X2, Y2 = registers.gens
    registers.add(X1 + X2)
    result = dcp.solve_multi_dcp(k, l, glv, registers)
    assert int(result[0]) in [
        744737971283560328960748838057,
        266501786421819480514269398178,
        192522633656806613321873142228,
    ]


def test_solve_multi_dcp_pari():
    k1, k2 = ZZ(2), ZZ(4)
    l = ZZ(3)
    glv = utils.GLVCurve()
    p = ZZ(1039)
    F = GF(p)
    glv.curve = EllipticCurve(F, [0, 6])
    P = glv.curve(136, 830)
    Q = 782 * P
    glv.order = ZZ(1033)
    glv.lam = ZZ(195)
    glv.beta = ZZ(140)
    assert (l * P)[0] + Q[0] - 54 == 0
    registers = utils.Registers()
    X1, Y1, X2, Y2 = registers.gens
    registers.add(X1 + X2 - 54)
    result = dcp.solve_glv_multi_dcp_pari([k1, k2], l, glv, registers)
    assert result[0] == P[0]


def test_solve_multi_dcp_pari2():
    k1, k2 = ZZ(2), ZZ(4)
    l = ZZ(-3)
    glv = utils.GLVCurve()
    p = ZZ(1039)
    F = GF(p)
    glv.curve = EllipticCurve(F, [0, 6])
    P = glv.curve(413, 459)
    Q = 782 * P
    glv.order = ZZ(1033)
    glv.lam = ZZ(195)
    glv.beta = ZZ(140)
    assert (l * P)[0] + Q[0] - 461 == 0
    registers = utils.Registers()
    X1, Y1, X2, Y2 = registers.gens
    registers.add(X1 + X2 - 461)
    result = dcp.solve_glv_multi_dcp_pari([k1, k2], l, glv, registers)
    assert result[0] == P[0]


def test_bin():
    assert msm.to_bin(2, 2) == [0, 1]
    assert msm.from_bin([0, 1, 1]) == 6
    assert 6 == msm.from_bin(msm.to_bin(6, 5))
    assert msm.zero_allign_bin(3, 3) == "0b011"


def test_twodim_recoding():
    l = 5
    kbin, kpbin = msm.to_bin(1, l), msm.to_bin(0, l)
    bbin, bpbin = msm.twodim_recoding(kbin, kpbin)
    bbin2, bpbin2 = msm.twodim_recoding2(kbin, kpbin)
    assert bbin == (l - 1) * [-1] + [1]
    assert bpbin == [0] * l
    assert bbin2 == (l - 1) * [-1] + [1]
    assert bpbin2 == [0] * l
    assert kbin, kpbin == msm.twodim_recoding_inverse(bbin, bpbin)

    kbin, kpbin = [1, 0, 0, 1, 1], [1, 0, 1, 1, 0]
    bbin, bpbin = msm.twodim_recoding(kbin, kpbin)
    bbin2, bpbin2 = msm.twodim_recoding2(kbin, kpbin)
    assert bbin == [-1, -1, 1, 1, 1]
    assert bpbin == [-1, -1, 0, 0, 1]
    assert bbin2 == [-1, -1, 1, 1, 1]
    assert bpbin2 == [-1, -1, 0, 0, 1]
    assert kbin, kpbin == msm.twodim_recoding_inverse(bbin, bpbin)


def test_two_dim_recoding_inverse_upper_bits():
    bbin = [-1, -1, 1, 1, 1][1:]
    bpbin = [-1, -1, 0, 0, 1][1:]
    l = 5
    assert msm.two_dim_recoding_inverse_upper_bits(l, bbin, bpbin) == (
        [0, 1, 1],
        [[0, 1, 0], [1, 1, 0]],
    )
    l = 16
    bbin = [1, 1, 1, -1, 1, 1]
    bpbin = [0, 0, 0, -1, 0, 1]
    assert msm.two_dim_recoding_inverse_upper_bits(l, bbin, bpbin) == (
        [1, 1, 1, 0, 1],
        [[1, 1, 0, 1, 0], [0, 0, 1, 1, 0]],
    )


def test_twodim_recoding_extended():
    l = 5
    bbin, bpbin = msm.twodim_recoding_extended(
        msm.to_bin(1, l), -1, msm.to_bin(0, l), 1
    )
    assert bbin == (l - 1) * [1] + [-1], bbin
    assert bpbin == [0] * l


def test_scalar_mul():
    E = EllipticCurve(GF(31), [0, 3])
    assert msm.scalar_mul(23, E(25, 2), lamb=36, beta=25, n=43) == E(14, 22)
    assert msm.scalar_mul_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43) == E(9, 9)
    assert msm.interleaving_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43, w=2) == E(
        9, 9
    )
    assert msm.interleaving_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43, w=3) == E(
        9, 9
    )
    assert msm.interleaving_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43, w=4) == E(
        9, 9
    )
    assert msm.interleaving_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43, w=5) == E(
        9, 9
    )

    p = ZZ(285443622588986402140519271645217346819)
    F = GF(p)
    curve = EllipticCurve(F, [0, 2])
    k0, k1 = ZZ(10082204432554003063), ZZ(14822413726714981050)
    lam, order = ZZ(253203832052255696853947885031684638103), ZZ(
        285443622588986402111823706003450522417
    )
    beta = ZZ(3783202341796351462745182158696532965)
    P = curve(
        57048526346403197948054049175117134356, 235690989991545003732641645973152703971
    )
    assert (
        msm.interleaving_positive(k0, k1, P, lamb=lam, beta=beta, n=order, w=2)
        == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_positive(k0, k1, P, lamb=lam, beta=beta, n=order, w=3)
        == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_positive(k0, k1, P, lamb=lam, beta=beta, n=order, w=4)
        == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_positive(k0, k1, P, lamb=lam, beta=beta, n=order, w=5)
        == (k0 + k1 * lam) * P
    )

    ok0, ok1 = 2 * (k0 // 2) - 1, 2 * (k1 // 2) - 1
    assert (
        msm.signed_shamir_scalar_mul_positive(
            ok0, ok1, P, lamb=lam, beta=beta, n=order, w=1
        )
        == (ok0 + ok1 * lam) * P
    )

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k0), ZZ(k1)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.registers = TEST_REGISTERS

    assert (
        msm.interleaving_oracle_positive(P, zvpparams, 0, w=2)[0] == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_oracle_positive(P, zvpparams, 0, w=3)[0] == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_oracle_positive(P, zvpparams, 0, w=4)[0] == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_oracle_positive(P, zvpparams, 0, w=5)[0] == (k0 + k1 * lam) * P
    )

    assert (
        msm.interleaving_easy_oracle_positive(P, zvpparams, w=2)[0]
        == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_easy_oracle_positive(P, zvpparams, w=3)[0]
        == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_easy_oracle_positive(P, zvpparams, w=4)[0]
        == (k0 + k1 * lam) * P
    )
    assert (
        msm.interleaving_easy_oracle_positive(P, zvpparams, w=5)[0]
        == (k0 + k1 * lam) * P
    )

    secrets.k0, secrets.k1 = ZZ(ok0), ZZ(ok1)
    assert (
        msm.signed_shamir_scalar_mul_oracle_positive(P, zvpparams, 0, w=1)[0]
        == (ok0 + ok1 * lam) * P
    )

    assert (
        msm.regular_interleaving_easy_oracle_positive(P, zvpparams, w=1)[0]
        == (ok0 + ok1 * lam) * P
    )

    assert (
        msm.regular_interleaving_easy_oracle_positive(P, zvpparams, w=2)[0]
        == (ok0 + ok1 * lam) * P
    )

    assert (
        msm.regular_interleaving_easy_oracle_positive(P, zvpparams, w=3)[0]
        == (ok0 + ok1 * lam) * P
    )


def test_shamir_scalar_mul():
    E = EllipticCurve(GF(31), [0, 3])
    assert msm.shamir_scalar_mul_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43) == E(
        9, 9
    )


def test_scalar_mul_oracle_positive():

    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(29), ZZ(26)
    lam, order = ZZ(195), ZZ(1033)
    beta = ZZ(140)
    P = curve(47, 353)
    Q = 143 * P
    Q2 = (6 + 6 * lam) * P

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order

    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.registers = TEST_REGISTERS

    assert Q2 == Q
    assert P[0] + Q[0] + 2 == 0
    assert msm.scalar_mul_oracle_positive(P, zvpparams, 2)[1]


def test_scalar_mul_oracle_positive2():
    p = ZZ(2143)
    F = GF(p)
    curve = EllipticCurve(F, [0, 5])
    k1, k2 = ZZ(0b1010101), ZZ(0b100010)
    lam, order = ZZ(1262), ZZ(2089)
    beta = ZZ(1793)
    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(4)
    zvpparams.registers = TEST_REGISTERS

    P = curve(1050, 236)
    assert msm.scalar_mul_oracle_positive(P, zvpparams, 4)[1]


def test_scalar_mul_oracle_positive3():
    p = ZZ(4294968187)
    F = GF(p)
    curve = EllipticCurve(F, [0, 3])
    k = ZZ(2613661231)
    k1, k2 = ZZ(15788), ZZ(2360)
    lam = ZZ(2168662147)
    beta = ZZ(4037544512)
    order = ZZ(4295070547)
    P = curve(912661074, 590983812)

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(0)
    zvpparams.registers = TEST_REGISTERS

    assert (k1 + k2 * lam) % order == k
    assert (k1 + k2 * lam) * P == k * P
    assert msm.scalar_mul(k, P, lam, beta, order) == k * P
    assert msm.scalar_mul_positive(k1, k2, P, lam, beta, order) == k * P
    assert msm.scalar_mul_oracle_positive(P, zvpparams, 0)[0] == k * P


def test_scalar_mul_oracle_positive4():
    p = ZZ(1099511628079)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(0b111001110011100110001), ZZ(0b111111011010100101101)
    lam = ZZ(533962724094)
    beta = ZZ(35271754617)
    order = ZZ(1099512457333)
    P = curve(961630273290, 234600338115)
    k = (k1 + k2 * lam) % order

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.registers = TEST_REGISTERS

    assert msm.scalar_mul(k, P, lam, beta, order) == k * P
    assert msm.scalar_mul_positive(k1, k2, P, lam, beta, order) == k * P
    assert msm.scalar_mul_oracle_positive(P, zvpparams, 0)[0] == k * P


def test_shamir_scalar_mul_oracle_positive():

    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(29), ZZ(26)
    lam, order = ZZ(195), ZZ(1033)
    beta = ZZ(140)
    P = curve(47, 353)
    Q = 143 * P
    Q2 = (6 + 6 * lam) * P

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.registers = TEST_REGISTERS

    assert Q2 == Q
    assert P[0] + Q[0] + 2 == 0
    assert msm.shamir_scalar_mul_oracle_positive(P, zvpparams, 2)[1]


def test_shamir_scalar_mul_oracle_positive3():
    p = ZZ(4294968187)
    F = GF(p)
    curve = EllipticCurve(F, [0, 3])
    k = ZZ(2613661231)
    k1, k2 = ZZ(15788), ZZ(2360)
    lam = ZZ(2168662147)
    beta = ZZ(4037544512)
    order = ZZ(4295070547)
    P = curve(912661074, 590983812)

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.registers = TEST_REGISTERS

    assert (k1 + k2 * lam) % order == k
    assert (k1 + k2 * lam) * P == k * P
    assert msm.shamir_scalar_mul_oracle_positive(P, zvpparams, 0)[0] == k * P


def test_shamir_scalar_mul_oracle_positive4():
    p = ZZ(1099511628079)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(0b111001110011100110001), ZZ(0b111111011010100101101)
    lam = ZZ(533962724094)
    beta = ZZ(35271754617)
    order = ZZ(1099512457333)
    P = curve(961630273290, 234600338115)
    k = (k1 + k2 * lam) % order

    glv = utils.GLVCurve()
    glv.curve = curve
    glv.lam = lam
    glv.beta = beta
    glv.order = order
    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(k1), ZZ(k2)

    zvpparams = ZVPparams(glv.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.registers = TEST_REGISTERS

    assert msm.shamir_scalar_mul_oracle_positive(P, zvpparams, 0)[0] == k * P


def test_zvp_sac_attack():
    glv = utils.GLVCurve()
    p = ZZ(2143)
    F = GF(p)
    glv.curve = EllipticCurve(F, [0, 5])
    glv.order = ZZ(2089)
    glv.lam = ZZ(1262)
    glv.beta = ZZ(1793)

    secrets = utils.GLVSecrets(glv)
    secrets.k0, secrets.k1 = ZZ(0b1010101), ZZ(0b100010)

    zvpparams = ZVPparams(p.nbits())
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(4)
    zvpparams.attack = zvpmsm.zvp_glv_sac

    zvpparams.registers = TEST_REGISTERS
    zvpparams.do_attack()

    l = zvpparams.secrets.half_bits + 1
    b0, b1 = msm.twodim_recoding(msm.to_bin(secrets.k0, l), msm.to_bin(secrets.k1, l))
    assert [b0[-4:], b1[-4:]] in zvpparams.results["scalars"]


def test_zvp_sac_attack2():
    secrets = utils.GLVSecrets()
    secrets.from_dict(
        {
            "k0": 12030649258185014393,
            "k1": 12631350741957173217,
            "k": 35962488788555180667629352222831099190,
            "curve": [0, 3],
            "field": 258902860607755093864506483129717834787,
            "beta": 74377596896322477575958861898445004046,
            "lambda": 153433049894344823205796386208105161084,
            "order": 258902860607755093877601888835855633873,
        }
    )
    zvpparams = ZVPparams(secrets.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(1)
    zvpparams.attack = zvpmsm.zvp_glv_sac

    zvpparams.registers = TEST_REGISTERS
    TEST_REGISTERS.add(X1 + X2)
    TEST_REGISTERS.add(X1 + X2 + 1)
    TEST_REGISTERS.add(X1 + X2 + 2)
    TEST_REGISTERS.add(X1**2 + X1 * X2 + X2**2)

    zvpparams.do_attack()


def test_zvp_sac_attack3():
    secrets = utils.GLVSecrets()
    secrets.from_dict(
        {
            "k0": 9354505673776504235,
            "k1": 16488824988860969087,
            "k": 29909316888294194379575008222298337520,
            "curve": [0, 5],
            "field": 219367660988620325471586898010713677223,
            "beta": 46913353786251974458108963163358176822,
            "lambda": 92096968189742899780809622067096840422,
            "order": 219367660988620325450527123821457422649,
        }
    )
    zvpparams = ZVPparams(secrets.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.attack = zvpmsm.zvp_glv_sac

    registers = utils.Registers()
    X1, Y1, X2, Y2, A, B = registers.all_gens
    registers.add_tuple(
        *(X1**3 * X2**3 + X1**3 * B + X2**3 * B + B**2 - 1, Y1 * Y2 + 1)
    )
    registers.add_tuple(
        *(
            -(X1**3) * X2**3 - X1**3 * B - X2**3 * B + 8 * B**2,
            Y1 * Y2 - 3 * B,
        )
    )
    registers.add_tuple(*(X1**2 * X2**2 + X1 * B + X2 * B, Y1 * X2 + X1 * Y2))
    zvpparams.registers = registers
    zvpparams.do_attack()


def test_zvp_sac_attack4():
    secrets = utils.GLVSecrets()
    secrets.from_dict(
        {
            "k0": 15043661883176416155,
            "k1": 10882924498039070729,
            "k": 89991237852242282884350750717970618338,
            "curve": [0, 3],
            "field": 305569605697985967340135829209822005883,
            "beta": 146834822803581668919535781180527474338,
            "lambda": 99344104049440347209412087658597197467,
            "order": 305569605697985967371993398978446362347,
        }
    )
    zvpparams = ZVPparams(secrets.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.attack = zvpmsm.zvp_glv_sac

    X1, Y1, X2, Y2, A, B = zvpparams.registers.all_gens
    zvpparams.registers.add_tuple(
        *(X1**2 * X2**2 + X1 * B + X2 * B, Y1 * X2 + X1 * Y2)
    )
    zvpparams.registers.add_tuple(
        *(
            X1**3 * X2**3 - X1**2 * X2**2 + X1**3 * B + X2**3 * B + B**2,
            X1 * X2 + Y1 * Y2,
        )
    )
    zvpparams.registers.add_tuple(
        *(
            -4 * X1**2 * X2**2 - X2**4 + 4 * X1 * B,
            X1**3 - 2 * X1**2 * X2 + X1 * X2**2 - Y1**2 + 2 * Y1 * Y2 - Y2**2,
        )
    )

    zvpparams.do_attack()


def test_glv_dcp_pari5():
    s0, s1 = 1, 0
    g0, g1 = -1, 0
    i = 63
    secrets = utils.GLVSecrets()
    secrets.from_dict(
        {
            "k0": 15043661883176416155,
            "k1": 10882924498039070729,
            "k": 89991237852242282884350750717970618338,
            "curve": [0, 3],
            "field": 305569605697985967340135829209822005883,
            "beta": 146834822803581668919535781180527474338,
            "lambda": 99344104049440347209412087658597197467,
            "order": 305569605697985967371993398978446362347,
        }
    )
    zvpparams = ZVPparams(secrets.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.attack = zvpmsm.zvp_glv_sac

    X1, Y1, X2, Y2, A, B = zvpparams.registers.all_gens
    zvpparams.registers.add_tuple(
        *(
            -4 * X1**2 * X2**2 - X2**4 + 4 * X1 * B,
            X1**3 - 2 * X1**2 * X2 + X1 * X2**2 - Y1**2 + 2 * Y1 * Y2 - Y2**2,
        )
    )

    dcp_scalars = 2 * s0 * g0 + 2 * s1 * g1 - 2 * s0 * g1, 2 * s1 * g0 - 2 * s0 * g1
    assert dcp_scalars == (-2, 0)

    P = dcp.solve_dcp(-2, zvpparams.glv, zvpparams.registers)


def test_Vpolynomial_degree():
    secrets = utils.GLVSecrets()
    secrets.from_dict(
        {
            "k0": 15043661883176416155,
            "k1": 10882924498039070729,
            "k": 89991237852242282884350750717970618338,
            "curve": [0, 3],
            "field": 305569605697985967340135829209822005883,
            "beta": 146834822803581668919535781180527474338,
            "lambda": 99344104049440347209412087658597197467,
            "order": 305569605697985967371993398978446362347,
        }
    )
    zvpparams = ZVPparams(secrets.bits)
    zvpparams.secrets = secrets
    zvpparams.target_bits = ZZ(2)
    zvpparams.attack = zvpmsm.zvp_glv_sac

    X1, Y1, X2, Y2, A, B = zvpparams.registers.all_gens
    zvpparams.registers.add_tuple(
        *(X1**3 * X2**3 + X1**3 * B + X2**3 * B + B**2 - 1, Y1 * Y2 + 1)
    )
    zvpparams.registers.add_tuple(
        *(
            -(X1**3) * X2**3 - X1**3 * B - X2**3 * B + 8 * B**2,
            Y1 * Y2 - 3 * B,
        )
    )
    zvpparams.registers.add_tuple(
        *(X1**2 * X2**2 + X1 * B + X2 * B, Y1 * X2 + X1 * Y2)
    )
    zvpparams.registers.add_tuple(
        *(
            X1**3 * X2**3 - X1**2 * X2**2 + X1**3 * B + X2**3 * B + B**2,
            X1 * X2 + Y1 * Y2,
        )
    )
    zvpparams.registers.add_tuple(
        *(
            -4 * X1**2 * X2**2 - X2**4 + 4 * X1 * B,
            X1**3 - 2 * X1**2 * X2 + X1 * X2**2 - Y1**2 + 2 * Y1 * Y2 - Y2**2,
        )
    )
    zvpparams.registers.add_tuple(
        *(
            -(X1**4) + 4 * X1**3 * X2 + 8 * X1 * B + 4 * X2 * B,
            2 * X1**3 - 3 * X1**2 * X2 + X2**3 - Y1**2 + 2 * Y1 * Y2 - Y2**2,
        )
    )

    Vpolynomials = [
        prepare_Vpolynomial(polynomial, zvpparams.registers, secrets.glv)
        for polynomial in zvpparams.registers.polynomials
    ]
    for V in Vpolynomials:
        print([f for f, _ in V.factor() if f in ZZ], V.degree(), len(V.monomials()))


def test_wnaf():
    # Example 3.34
    k = 1122334455
    assert msm.wnaf(k, 2) == [
        -1,
        0,
        0,
        -1,
        0,
        0,
        0,
        0,
        -1,
        0,
        0,
        -1,
        0,
        0,
        0,
        -1,
        0,
        -1,
        0,
        1,
        0,
        -1,
        0,
        0,
        -1,
        0,
        1,
        0,
        0,
        0,
        1,
    ]
    assert msm.wnaf(k, 3) == [
        -1,
        0,
        0,
        -1,
        0,
        0,
        0,
        0,
        -1,
        0,
        0,
        -1,
        0,
        0,
        0,
        3,
        0,
        0,
        1,
        0,
        0,
        -1,
        0,
        0,
        3,
        0,
        0,
        0,
        0,
        0,
        1,
    ]
    assert msm.wnaf(k, 4) == [
        7,
        0,
        0,
        0,
        -1,
        0,
        0,
        0,
        7,
        0,
        0,
        0,
        7,
        0,
        0,
        0,
        5,
        0,
        0,
        0,
        0,
        7,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        1,
    ]
    assert msm.wnaf(k, 5) == [
        -9,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        0,
        0,
        11,
        0,
        0,
        0,
        0,
        0,
        -9,
        0,
        0,
        0,
        0,
        -15,
        0,
        0,
        0,
        0,
        1,
    ]

    assert msm.from_wnaf(msm.wnaf(k, 2)) == k
    assert msm.from_wnaf(msm.wnaf(k, 3)) == k
    assert msm.from_wnaf(msm.wnaf(k, 4)) == k
    assert msm.from_wnaf(msm.wnaf(k, 5)) == k


def test_regular_wnaf():
    # Example 3.34
    k = 1122334455
    # k = 13
    assert msm.from_regular_wnaf(msm.regular_wnaf(k, 1), 1) == k
    assert msm.from_regular_wnaf(msm.regular_wnaf(k, 2), 2) == k
    assert msm.from_regular_wnaf(msm.regular_wnaf(k, 3), 3) == k
    assert msm.from_regular_wnaf(msm.regular_wnaf(k, 4), 4) == k
    assert msm.from_regular_wnaf(msm.regular_wnaf(k, 5), 5) == k


def test_bit_length():
    for k in range(1, 10000, 2):
        assert len(msm.regular_wnaf(k, 3)) == ceil(k.bit_length() / 3), k
    for _ in range(100):
        k = randint(2**25, 2**30)
        k = 2 * k + 1
        assert len(msm.regular_wnaf(k, 3)) == ceil(k.bit_length() / 3), k
        assert len(msm.regular_wnaf(k, 2)) == ceil(k.bit_length() / 2), k
        assert len(msm.regular_wnaf(k, 4)) == ceil(k.bit_length() / 4), k
        assert len(msm.regular_wnaf(k, 5)) == ceil(k.bit_length() / 5), k
