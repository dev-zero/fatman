
export interface TestResult {
  id: number;
  ctime: string;
  data: { [id: string]: Object; };
  test: { id: number; name: string; };
  method: { id: number; pseudopotential: string };
}

export function toTestResult(r: any) : TestResult {
  return <TestResult>({
    id: r.id,
    ctime: r.ctime,
    data: r.data,
    test: r.test,
    method: r.method
  });
}
