import classNames from 'classnames';
import { ReactNode } from 'react';
import { ResetType } from 'solarxr-protocol';
import { useLayout } from '../hooks/layout';
import { BVHButton } from './BVHButton';
import { ResetButton } from './home/ResetButton';
import { Navbar } from './Navbar';
import { TopBar } from './TopBar';
import { OverlayWidget } from './widgets/OverlayWidget';

export function MainLayoutRoute({
  children,
  background = true,
  widgets = true,
}: {
  children: ReactNode;
  background?: boolean;
  widgets?: boolean;
}) {
  const { layoutHeight, ref } = useLayout<HTMLDivElement>();
  const { layoutWidth, ref: refw } = useLayout<HTMLDivElement>();

  return (
    <>
      <TopBar></TopBar>
      <div ref={ref} className="flex-grow" style={{ height: layoutHeight }}>
        <div className="flex h-full pb-3 flex-col xs:flex-row gap-2 xs:gap-0">
          <Navbar></Navbar>
          <div
            className="flex gap-2 pr-3 w-full flex-col-reverse xs:flex-row px-2 xs:px-0"
            ref={refw}
          >
            <div
              className={classNames(
                'flex flex-col rounded-xl w-full overflow-hidden',
                background && 'bg-background-70'
              )}
            >
              {children}
            </div>
            {widgets && (
              <div className="flex flex-row xs:flex-col px-2 xs:min-w-[274px] xs:w-[274px] gap-2 pt-2 pb-2 xs:pb-0 rounded-xl overflow-y-auto bg-background-70">
                <div className="grid grid-cols-2 gap-2 w-full">
                  <ResetButton type={ResetType.Quick}></ResetButton>
                  <ResetButton type={ResetType.Full}></ResetButton>
                  <ResetButton type={ResetType.Mounting}></ResetButton>
                  <BVHButton></BVHButton>
                </div>
                <div className="w-full">
                  <OverlayWidget></OverlayWidget>
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </>
  );
}
