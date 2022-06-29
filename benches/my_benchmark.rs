use castepbin_rust::dos::dos_util::calculate_total_dos_from_bands;
use criterion::{criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(c: &mut Criterion) {
    let band_file = "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/castepbin-rust/Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.bands";
    c.bench_function("total dos", |b| {
        b.iter(|| calculate_total_dos_from_bands(band_file, None, None))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
